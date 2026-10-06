"""Build the ClinVar variant list for the clinical-label part of the study.

Reads one ClinVar release (the outcome release), the ClinVar files that define what was already labelled
(the developers' 2025-06-15 VCF and the 2025-07 and 2024-07 variant_summary archives) and the derived
exclusion files, and writes one row per pathogenic/likely pathogenic or benign/likely benign variant of
the outcome release:

    <out>/clinvar_variants.tsv.gz   one row per variant, with its region and set membership
    <out>/clinvar_variants_summary.json   row flow, exclusions and counts per set, region and chromosome group

A variant is in the new-label set when it has no earlier pathogenic/benign classification in any history
file (by allele, AlleleID or VariationID through ClinVar's replacement history), is not at a position of
the developers' benchmark data, has no submission from the organisation that supplied their GREGoR cases,
and its identity is resolved. Every variant is in the all-labels set.

Usage:
    python src/variant_sets.py --release <variant_summary> --e-dev <vcf> --t0 <variant_summary>
        --e-hist <variant_summary> --submission-summary <gz> [--out data/derived/variants]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from collections import Counter
from pathlib import Path

import clinvar_index as idx
import clinvar_records as cr
import hgvs_strata as hs
import liftover
from paths import DERIVED, EXPOSURE, sha256_file

FILES = {
    "positions": EXPOSURE / "development_positions.tsv",
    "submitters": EXPOSURE / "development_submitters.txt",
    "replacements": EXPOSURE / "clinvar_replacements.tsv",
    "mane": EXPOSURE / "mane_select.tsv",
    "chain": EXPOSURE / "hg19ToHg38.over.chain.gz",
}
SUBMITTER_SOURCE = "development submitters"
COLUMNS = ("variation_id", "allele_id", "chrom", "pos", "ref", "alt", "gene_id", "label", "stars",
           "review_status", "name", "chromosome_group", "variant_type", "family", "region", "intron_offset",
           "mane_select", "new_label", "exclusion", "development_material", "labelled_before_in")


def read_submitter_vids(summary: Path, names_file: Path, flow: Counter) -> set:
    """VariationIDs with a submission from an organisation listed in `names_file` (TSV: organization_id,
    submitter_name). Stops if a listed name matches no submission, so a changed spelling cannot silently
    empty the exclusion."""
    lines = names_file.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0].split("\t") != ["organization_id", "submitter_name"]:
        raise idx.InputFormatError(f"{names_file.name}: header must be organization_id, submitter_name")
    names = {}
    for line in lines[1:]:
        if line.strip():
            org, name = line.split("\t")
            names[name.strip()] = org
            names[name.strip().replace("; ", ", ")] = org
    vids, matched = set(), Counter()
    for vid, submitter in idx.read_submitters(summary, flow):
        if submitter in names and vid.isdigit():
            vids.add(vid)
            matched[names[submitter]] += 1
    for org in sorted(set(names.values())):
        if not matched[org]:
            raise idx.InputFormatError(f"organisation {org} matches no submission")
        flow[f"submission_summary: submissions from organisation {org}"] = matched[org]
    return vids


def read_mane(path: Path) -> set:
    """MANE Select RefSeq transcripts (first column of the derived table)."""
    lines = path.read_text(encoding="utf-8").splitlines()
    return {line.split("\t")[0] for line in lines[1:] if line.strip()}


def name_transcript(name: str) -> str | None:
    head = (name or "").split("(", 1)[0].split(":", 1)[0]
    return head if head.startswith(("NM_", "NR_")) and "." in head else None


def exclusion(c: idx.Candidate) -> str:
    """The first reason a variant is not in the new-label set, or ''."""
    if c.unresolved:
        return c.unresolved
    if c.development:
        return "developers' benchmark or submitter"
    if c.seen:
        return "labelled before"
    return ""


def row(c: idx.Candidate, mane: set) -> dict:
    rec = c.record
    chrom, pos, ref, alt = rec.key.split(":")
    assignment = hs.assign(rec.name)
    reason = exclusion(c)
    return {
        "variation_id": rec.variation_id, "allele_id": rec.allele_id, "chrom": chrom, "pos": pos,
        "ref": ref, "alt": alt, "gene_id": rec.gene, "label": c.group, "stars": rec.stars,
        "review_status": rec.review_status, "name": rec.name,
        "chromosome_group": cr.chromosome_group(rec.chrom), "variant_type": cr.variant_type(rec.key),
        "family": assignment.family, "region": assignment.stratum,
        "intron_offset": "" if assignment.signed_offset is None else assignment.signed_offset,
        "mane_select": "yes" if name_transcript(rec.name) in mane else "no",
        "new_label": "no" if reason else "yes", "exclusion": reason,
        "development_material": c.development or "",
        "labelled_before_in": ",".join(sorted(c.seen)),
    }


def build(release: Path, exposure_paths: dict, submission_summary: Path, out: Path) -> dict:
    flow = Counter()
    dev_vids = read_submitter_vids(submission_summary, FILES["submitters"], flow)
    aliases = idx.read_replacements(FILES["replacements"])
    mane = read_mane(FILES["mane"])
    lift = liftover.load(FILES["chain"])
    universe = idx.T1Universe(canon=lambda v: aliases.get(v, v) if v is not None else None)
    candidates = idx.assemble_t1(idx.read_variant_summary(release, "T1", flow, fail_closed=False), flow, universe)
    index = idx.ExposureIndex(candidates, universe, aliases)
    for source in idx.BASE_EXPOSURE:
        reader = idx.read_vcf if source == "E-dev" else idx.read_variant_summary
        for rec in reader(exposure_paths[source], source, flow):
            index.add(rec, flow)
    positions = idx.read_positions(FILES["positions"])
    development = idx.mark_development(candidates, positions, {SUBMITTER_SOURCE: dev_vids}, index.canon)
    del positions
    blocked = index.block_unmatched(lift)

    out.mkdir(parents=True, exist_ok=True)
    counts = Counter()
    with gzip.open(out / "clinvar_variants.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for vid in sorted(candidates, key=lambda v: int(v) if v.isdigit() else -1):
            r = row(candidates[vid], mane)
            writer.writerow(r)
            stars = r["stars"] if isinstance(r["stars"], int) else "unknown"
            for label_set in ("all labels", "new labels") if r["new_label"] == "yes" else ("all labels",):
                counts[(label_set, r["family"], r["region"], r["chromosome_group"], r["variant_type"],
                        stars, r["label"])] += 1
    summary = {
        "inputs": {str(p): sha256_file(p) for p in (release, *exposure_paths.values(), submission_summary,
                                                    *FILES.values())},
        "flow": dict(sorted(flow.items())),
        "development_material": dict(development),
        "unknown_record_blocking": blocked,
        "counts": [{"set": k[0], "family": k[1], "region": k[2], "chromosome_group": k[3],
                    "variant_type": k[4], "stars": k[5], "label": k[6], "n": v}
                   for k, v in sorted(counts.items(), key=lambda kv: str(kv[0]))],
    }
    (out / "clinvar_variants_summary.json").write_text(json.dumps(summary, indent=1, default=str) + "\n",
                                                       encoding="utf-8", newline="\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    for name in ("--release", "--e-dev", "--t0", "--e-hist", "--submission-summary"):
        ap.add_argument(name, required=True, type=Path)
    ap.add_argument("--out", type=Path, default=DERIVED / "variants")
    a = ap.parse_args(argv)
    build(a.release, {"E-dev": a.e_dev, "T0": a.t0, "E-hist": a.e_hist}, a.submission_summary, a.out)
    for name in ("clinvar_variants.tsv.gz", "clinvar_variants_summary.json"):
        print(f"{a.out / name}\tsha256 {sha256_file(a.out / name)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
