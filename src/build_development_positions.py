"""Development positions (see data/development_sources.md).

Writes data/derived/exposure/development_positions.tsv (chrom, start, end, source; GRCh38, 1-based,
merged per source) from the recorded preparation files:
  - saturation genome editing: MaveDB positions per study (src/mavedb_sources.py), or the study's
    gene region (RefSeq GFF named by MANE v1.5, padded by 1,000 bases, per chromosome) when
    the source table fixes the gene region, a file is missing, or any scored variant has no position;
  - TraitGym Mendelian and Complex: every row of the seven recorded test files (any label);
  - ncVarDB: every row of the four recorded files (hg38 only).
Then writes one evidence file per source and development_sources.tsv. No score or label is kept.
Outputs are never overwritten.

Usage: python src/build_development_positions.py
"""
from __future__ import annotations

import gzip
import json
import re
import sys
from collections import Counter
from pathlib import Path

import clinvar_records as cr
import mavedb_sources as ms
from build_exposure_files import write_sources
from paths import ROOT, sha256_file
from paths import RAW as RAW_DIR
from clinvar_index import InputFormatError

RAW = RAW_DIR
EXPOSURE = ROOT / "data" / "derived" / "exposure"
GFF = RAW / "GCF_000001405.40_GRCh38.p14_genomic.gff.gz"
TRAITGYM = {"TraitGym Mendelian": ("mendelian_traits_all", "mendelian_traits_matched_9",
                                   "mendelian_traits_v21_matched_9", "mendelian_traits_v22_matched_9"),
            "TraitGym Complex": ("complex_traits_all", "complex_traits_matched_9", "complex_traits_v22_matched_9")}
NCVAR_TSV = ("ncVar_pathogenic.tsv", "ncVar_benign.tsv")
NCVAR_VCF = ("ncVar_pathogenic.vcf", "ncVar_benign.vcf")
PAD = 1000
_ATTR = re.compile(r"([^;=]+)=([^;]*)")


def _attributes(text: str) -> dict:
    return {k: v for k, v in _ATTR.findall(text)}


def _gff_lines(path: Path):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for line in handle:
            if not line.startswith("#"):
                parts = line.rstrip("\r\n").split("\t")
                if len(parts) == 9:
                    yield parts


def chromosome_map(gff: Path = GFF) -> dict:
    """RefSeq accession -> chromosome name for the primary assembly chromosomes."""
    out = {}
    for seqid, _, kind, _, _, _, _, _, attrs in _gff_lines(gff):
        if kind == "region":
            a = _attributes(attrs)
            if a.get("genome") in ("chromosome", "mitochondrion") and a.get("chromosome"):
                out[seqid] = cr.canonical_chrom(a["chromosome"])
    if not out:
        raise InputFormatError(f"{gff.name}: no chromosome region")
    return out


def gene_regions(gff: Path, symbols, chroms: dict | None = None) -> dict:
    """symbol -> GeneID and one padded span per chromosome of its `gene` features (primary assembly)."""
    chroms = chroms or chromosome_map(gff)
    wanted, found = set(symbols), {}
    for seqid, _, kind, start, end, _, _, _, attrs in _gff_lines(gff):
        if kind != "gene" or seqid not in chroms:
            continue
        a = _attributes(attrs)
        if a.get("gene") not in wanted:
            continue
        ids = [x.split(":", 1)[1] for x in a.get("Dbxref", "").split(",") if x.startswith("GeneID:")]
        if len(ids) != 1:
            raise InputFormatError(f"gene feature for {a['gene']} without one GeneID")
        spans = found.setdefault(a["gene"], {}).setdefault(ids[0], {})
        lo, hi = spans.get(chroms[seqid], (int(start), int(end)))
        spans[chroms[seqid]] = (min(lo, int(start)), max(hi, int(end)))
    out = {}
    for symbol in symbols:
        genes = found.get(symbol, {})
        if len(genes) != 1:
            raise InputFormatError(f"{symbol}: {len(genes)} GeneIDs on the primary assembly (one required)")
        (gene_id, spans), = genes.items()
        out[symbol] = {"gene_id": gene_id, "intervals": sorted((c, max(1, s - PAD), e + PAD) for c, (s, e) in spans.items())}
    return out


def span(chrom: str, pos, ref: str, alt: str) -> tuple:
    """Full span of a variant (source table: any allele, full span for indels).
    VCF style: the reference bases; for an insertion (an alternate allele longer than the reference)
    also the base after it, so that both flanks are covered. ncVarDB '-' style insertion: the bases
    before, at and after `pos`, which covers both flanks whichever base `pos` names."""
    p = int(pos)
    chrom = cr.canonical_chrom(str(chrom))
    if ref == "-":
        return chrom, p - 1, p + 1
    end = p + len(ref) - 1
    if any(a not in ("-", "", ".") and len(a) > len(ref) for a in (alt or "").split(",")):
        end = p + len(ref)
    return chrom, p, end


def traitgym_intervals(path: Path) -> list:
    """Span of every row (positive or control)."""
    import pyarrow.parquet as pq
    table = pq.read_table(path, columns=["chrom", "pos", "ref", "alt"])
    out = []
    for chrom, pos, ref, alt in zip(*(table.column(c).to_pylist() for c in ("chrom", "pos", "ref", "alt"))):
        if chrom is None or not isinstance(pos, int) or not ref or alt is None:
            raise InputFormatError(f"{path}: row without chrom, pos, ref or alt")
        out.append(span(chrom, pos, ref, alt))
    return out


def ncvardb_intervals(tsvs, vcfs) -> list:
    """Every pathogenic and benign variant of ncVarDB (hg38 only). The TSV files declare the genome of
    every row; a VCF row is accepted only if its ncVarDB identifier (INFO NC) names a TSV row at the
    same chromosome and position, so that its assembly is established, and a declared VCF reference
    must be hg38 or GRCh38."""
    out, located = set(), {}
    for path in tsvs:
        lines = path.read_text(encoding="utf-8").splitlines()
        header = lines[0].split("\t")
        col = {f: header.index(f) for f in ("ID", "genome", "chr", "pos", "ref", "alt")}
        for line in lines[1:]:
            v = line.split("\t")
            if v[col["genome"]] != "hg38":
                raise InputFormatError(f"{path.name}: a row is not on hg38")
            here = (cr.canonical_chrom(v[col["chr"]]), int(v[col["pos"]]))
            for ident in (i.strip() for i in v[col["ID"]].split(",") if i.strip()):
                if located.setdefault(ident, here) != here:
                    raise InputFormatError(f"{path.name}: identifier {ident} at two positions")
            out.add(span(v[col["chr"]], v[col["pos"]], v[col["ref"]], v[col["alt"]]))
    for path in vcfs:
        for line in path.read_text(encoding="utf-8").splitlines():
            if line.lower().startswith("##reference=") and not any(
                    a in line.lower() for a in ("hg38", "grch38")):
                raise InputFormatError(f"{path.name}: declared reference is not hg38")
            if not line or line.startswith("#"):
                continue
            v = line.split("\t")
            info = dict(x.split("=", 1) for x in v[7].split(";") if "=" in x) if len(v) > 7 else {}
            ids = [i for i in info.get("NC", "").split(",") if i]
            here = (cr.canonical_chrom(v[0]), int(v[1]))
            if not ids or any(located.get(i) != here for i in ids):
                raise InputFormatError(f"{path.name}: row at {v[0]}:{v[1]} does not match its hg38 ncVarDB row")
            out.add(span(v[0], v[1], v[3], v[4]))
    return sorted(out)


def sge_choice(score_sets, region: dict, missing=()) -> tuple:
    """('positions', intervals) if no file is missing and every score set is complete and positioned,
    else ('gene region', ...)."""
    if missing or not score_sets or any(s["unpositioned"] or not s["scored"] for s in score_sets):
        return "gene region", list(region["intervals"])
    return "positions", sorted({iv for s in score_sets for iv in s["intervals"]})


def merged_rows(by_source: dict) -> list:
    rows = []
    for source, intervals in by_source.items():
        merged = []
        for chrom, start, end in sorted((cr.canonical_chrom(c), s, e) for c, s, e in intervals):
            if merged and merged[-1][0] == chrom and start <= merged[-1][2] + 1:
                merged[-1][2] = max(merged[-1][2], end)
            else:
                merged.append([chrom, start, end])
        rows += [(c, s, e, source) for c, s, e in merged]
    return rows


def write_positions(rows: list, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "x", encoding="utf-8", newline="\n") as handle:
        handle.write("chrom\tstart\tend\tsource\n")
        handle.writelines(f"{c}\t{s}\t{e}\t{src}\n" for c, s, e, src in rows)


def _write_evidence(name: str, content: dict) -> str:
    rel = f"evidence/{name}.json"
    path = EXPOSURE / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "x", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(content, indent=1, sort_keys=True) + "\n")
    return rel


def _inputs(paths) -> dict:
    return {str(Path(p).relative_to(ROOT)).replace("\\", "/"): sha256_file(Path(p)) for p in paths}


def build() -> dict:
    chroms = chromosome_map(GFF)
    regions = gene_regions(GFF, sorted({s["gene"] for s in ms.STUDIES}), chroms)
    sge, sge_files, studies = [], [GFF], []
    for study in ms.STUDIES:
        score_sets, missing, urns = None, [], []
        if study["experiments"] is not None:
            score_sets = []
            try:
                urns = ms.study_score_sets(study)
            except ms.MissingSource as gone:          # a record never obtained: gene region
                urns, missing = [], [str(gone)]
            for urn in urns:
                files = [RAW / ms.local_name(urn, k) for k in ("scores", "mapped")]
                missing += [ms.local_name(urn, k) for k, f in zip(("scores", "mapped"), files) if not f.exists()]
                if all(f.exists() for f in files):
                    sge_files += files
                    score_sets.append(ms.score_set_positions(*files, chroms))
        choice, intervals = sge_choice(score_sets, regions[study["gene"]], missing)
        sge += intervals
        studies.append({"study": study["study"], "gene": study["gene"], "gene_id": regions[study["gene"]]["gene_id"],
                        "used": choice, "score_sets": len(score_sets or []),
                        "unpositioned": sum(s["unpositioned"] for s in score_sets or []), "missing": missing,
                        "score_set_urns": urns if study["experiments"] is not None else [],
                        "post_mapped_forms": dict(sorted(sum((Counter(s["forms"]) for s in score_sets or []),
                                                             Counter()).items()))})
    by_source = {"saturation genome editing": sge}
    tg_files = {}
    for source, dirs in TRAITGYM.items():
        tg_files[source] = [RAW / "traitgym" / d / "test.parquet" for d in dirs]
        by_source[source] = [iv for p in tg_files[source] for iv in traitgym_intervals(p)]
    nc_tsv, nc_vcf = [RAW / "ncvardb" / n for n in NCVAR_TSV], [RAW / "ncvardb" / n for n in NCVAR_VCF]
    by_source["ncVarDB"] = ncvardb_intervals(nc_tsv, nc_vcf)
    rows = merged_rows(by_source)
    write_positions(rows, EXPOSURE / "development_positions.tsv")
    counts = {s: sum(1 for r in rows if r[3] == s) for s in by_source}
    evidence = {
        "saturation genome editing": _write_evidence("saturation_genome_editing", {
            "studies": studies, "inputs": _inputs(sge_files), "merged_intervals": counts["saturation genome editing"]}),
        "TraitGym Mendelian": _write_evidence("traitgym_mendelian", {
            "inputs": _inputs(tg_files["TraitGym Mendelian"]), "merged_intervals": counts["TraitGym Mendelian"]}),
        "TraitGym Complex": _write_evidence("traitgym_complex", {
            "inputs": _inputs(tg_files["TraitGym Complex"]), "merged_intervals": counts["TraitGym Complex"]}),
        "ncVarDB": _write_evidence("ncvardb", {"inputs": _inputs(nc_tsv + nc_vcf), "merged_intervals": counts["ncVarDB"]}),
        "GREGoR": _write_evidence("gregor", {
            "rule": "submitter rule", "submitter_file": "development_submitters.txt",
            "inputs": _inputs([RAW / "organization_summary.txt"])}),
        "complex-trait benchmark": _write_evidence("complex_trait_benchmark", {
            "rule": "not excluded (evaluated only on trait association; variant list not published)"}),
    }
    write_sources(evidence, EXPOSURE / "development_sources.tsv")
    return {"merged_intervals": counts, "studies": studies}


if __name__ == "__main__":
    print(json.dumps(build(), indent=1))
    sys.exit(0)
