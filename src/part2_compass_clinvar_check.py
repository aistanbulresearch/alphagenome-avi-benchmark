"""Provenance check of the COMPASS random-variant proxy (labels only, no scores).

The authors' 'clinvar_single' label mixes their ClinVar variants and random SNVs; the planned rule
(src/part2_data.py harmonise_compass) sets aside the alleles present in the 2025-06-15 ClinVar VCF and keeps the
rest as random. This check counts how many COMPASS main-set SNVs of each subset appear in the archived ClinVar
variant summaries of July 2024 and July 2025 (around the library's design) and of October 2026. A random-labelled
allele present in an archive from before the design would be a ClinVar variant that the proxy kept as random.
The same counts are also given for the alleles that enter the analysis (after the Part 2 exclusions of Part 1 new
labels and the developers' benchmark data), with the October 2026 ClinVar classification of the analysed random
alleles found there (they are not Part 1 new labels, so none should be P/LP or B/LB).

Output: results/part2/compass_clinvar_check.json (counts only).

Usage:
    python src/part2_compass_clinvar_check.py
"""
from __future__ import annotations

import csv
import gzip
import json
import sys
from pathlib import Path

from paths import DERIVED, RAW, ROOT

TABLE = DERIVED / "part2" / "assay_variants.tsv.gz"
RELEASES = ["variant_summary_2024-07.txt.gz", "variant_summary_2025-07.txt.gz", "variant_summary_2026-10.txt.gz"]
OUT = ROOT / "results" / "part2" / "compass_clinvar_check.json"


def norm(chrom: str) -> str:
    c = chrom[3:] if chrom.lower().startswith("chr") else chrom
    return "MT" if c == "M" else c


def compass_keys(path: Path = TABLE) -> dict:
    """{(chrom, pos, ref, alt): subset} of the COMPASS main-set SNVs."""
    keys = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if r["dataset"] == "compass" and r["use"] == "main" and r["variant_type"] == "SNV":
                keys[(norm(r["chrom"]), int(r["pos"]), r["ref"], r["alt"])] = r["subset"]
    return keys


def analysed_keys(keys: dict, part1: set, development: set) -> dict:
    """The keys whose allele ('chrom:pos:ref:alt') is neither a Part 1 new label nor in the developers' benchmark data."""
    return {k: v for k, v in keys.items()
            if f"{k[0]}:{k[1]}:{k[2]}:{k[3]}" not in part1 and f"{k[0]}:{k[1]}:{k[2]}:{k[3]}" not in development}


def present_by_subset(keys: dict, summary: Path) -> dict:
    """Number of keys of each subset found as GRCh38 alleles in a ClinVar variant_summary file."""
    hit = {}
    with gzip.open(summary, "rt", encoding="utf-8", errors="replace") as h:
        header = h.readline().lstrip("#").rstrip("\n").split("\t")
        ix = [header.index(k) for k in ("Assembly", "Chromosome", "PositionVCF", "ReferenceAlleleVCF",
                                        "AlternateAlleleVCF")]
        for line in h:
            f = line.rstrip("\n").split("\t")
            if f[ix[0]] != "GRCh38":
                continue
            try:
                k = (norm(f[ix[1]]), int(f[ix[2]]), f[ix[3]], f[ix[4]])
            except ValueError:
                continue
            if k in keys:
                hit[k] = keys[k]
    counts = {}
    for s in hit.values():
        counts[s] = counts.get(s, 0) + 1
    return counts


def significance_of(keys: dict, summary: Path, subset: str) -> dict:
    """ClinicalSignificance of the keys of one subset found in a variant_summary file, counted once per allele.

    Every row of an allele is read; an allele whose rows carry different classifications is counted under all of
    them joined with ' | ', so a later row cannot hide a P/LP or B/LB classification."""
    by_allele = {}
    with gzip.open(summary, "rt", encoding="utf-8", errors="replace") as h:
        header = h.readline().lstrip("#").rstrip("\n").split("\t")
        ix = [header.index(k) for k in ("Assembly", "Chromosome", "PositionVCF", "ReferenceAlleleVCF",
                                        "AlternateAlleleVCF", "ClinicalSignificance")]
        for line in h:
            f = line.rstrip("\n").split("\t")
            if f[ix[0]] != "GRCh38":
                continue
            try:
                k = (norm(f[ix[1]]), int(f[ix[2]]), f[ix[3]], f[ix[4]])
            except ValueError:
                continue
            if keys.get(k) == subset:
                by_allele.setdefault(k, set()).add(f[ix[5]])
    counts = {}
    for classes in by_allele.values():
        label = " | ".join(sorted(classes))
        counts[label] = counts.get(label, 0) + 1
    return counts


def _totals(keys: dict) -> dict:
    out = {}
    for s in keys.values():
        out[s] = out.get(s, 0) + 1
    return out


def main() -> int:
    import pandas as pd

    import part2_analysis as p2

    keys = compass_keys()
    frame = pd.DataFrame([{"allele": f"{c}:{p}:{r}:{a}", "chrom": c, "pos": str(p)} for c, p, r, a in keys])
    analysed = analysed_keys(keys, p2.part1_new_labels(), p2.development_hits(frame))
    out = {"compass_main_snvs": _totals(keys), "present_in": {r: present_by_subset(keys, RAW / r) for r in RELEASES},
           "analysed_snvs": _totals(analysed),
           "analysed_present_in": {r: present_by_subset(analysed, RAW / r) for r in RELEASES},
           "analysed_random_significance_2026_10": significance_of(analysed, RAW / RELEASES[-1], "random_snv")}
    OUT.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
