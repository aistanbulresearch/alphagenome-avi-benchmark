"""Per-variant scores of AVI, SpliceAI and CADD for the variants of the ClinVar variant list.

Each source is a bgzip file with a tabix index; only the blocks holding the requested positions are
read (src/gnomad_remote.py): AVI from the local copy of the AlphaGenome download, SpliceAI and CADD
from their public servers by byte ranges.

    AVI       AVI SNV scores (AlphaGenome downloads), raw score and Phred-scaled score
    SpliceAI  Ensembl MANE v1.4 precomputed file, masked, distance 50: maximum of DS_AG, DS_AL, DS_DG,
              DS_DL over the gene entries of the allele
    CADD      v1.7 whole-genome SNV file: RawScore and PHRED

Output: data/derived/scores/<tool>.tsv.gz (variation_id, chrom, pos, ref, alt, score columns) and
<tool>_source.json. SpliceAI and CADD scores are not redistributed; data/derived/scores/ is not in git.

Usage:
    python src/tool_scores.py --tool {avi,spliceai,cadd} [--avi-file <path to .bgz>] [--threads 24]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sys
import time
import urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Callable

import gnomad_frequencies as gf
import gnomad_remote as gr
from paths import DERIVED

OUT = DERIVED / "scores"
SPLICEAI = ("https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_plugins/"
            "spliceai_scores.masked.snv.ensembl_mane_v1.4.grch38.vcf.gz")
CADD = "https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz"


def spliceai_max(info: str, alt: str) -> str:
    """Maximum delta score over the SpliceAI entries of `alt` in a VCF INFO field, or ''."""
    best = None
    for item in info.split(";"):
        if not item.startswith("SpliceAI="):
            continue
        for entry in item[len("SpliceAI="):].split(","):
            f = entry.split("|")
            if len(f) >= 6 and f[0] == alt:
                values = [float(v) for v in f[2:6] if v not in ("", ".")]
                if values:
                    best = max(values) if best is None else max(best, *values)
    return "" if best is None else f"{best:g}"


def parse_spliceai(parts: list, ref: str, alt: str) -> dict | None:
    if parts[3] != ref or alt not in parts[4].split(","):
        return None
    score = spliceai_max(parts[7], alt)
    return {"spliceai_max": score} if score != "" else None


def parse_cadd(parts: list, ref: str, alt: str) -> dict | None:
    """CADD TSV: Chrom, Pos, Ref, Alt, RawScore, PHRED."""
    if parts[2] != ref or parts[3] != alt:
        return None
    return {"cadd_raw": parts[4], "cadd_phred": parts[5].strip()}


def parse_avi(columns: list) -> Callable:
    """AVI TSV with a header line naming its columns; the allele columns and the two scores are found
    by name in `columns` (the header of the file)."""
    lower = [c.lower().lstrip("#") for c in columns]

    def find(*names):
        for n in names:
            if n in lower:
                return lower.index(n)
        raise ValueError(f"AVI header lacks one of {names}: {columns}")
    i_ref, i_alt = find("ref", "reference"), find("alt", "alternate")
    i_raw, i_phred = find("avi", "avi_score", "score", "raw"), find("avi_phred", "phred", "avi_score_phred")

    def parse(parts: list, ref: str, alt: str) -> dict | None:
        fields = parts[:8] + (parts[8].split("\t") if len(parts) > 8 else [])
        if fields[i_ref] != ref or fields[i_alt] != alt:
            return None
        return {"avi_raw": fields[i_raw], "avi_phred": fields[i_phred].strip()}
    return parse


def local_fetch(path: Path) -> Callable:
    def fetch(_url, start, end):
        with open(path, "rb") as h:
            h.seek(start)
            return h.read(end - start)
    return fetch


def lookup_all(variants: list, url: str, index: gr.TabixIndex, contig_of: Callable, parse: Callable,
               fetch: Callable, threads: int) -> list:
    by_pos = defaultdict(list)
    for v in variants:
        by_pos[(v["chrom"], int(v["pos"]))].append(v)

    def one(key):
        chrom, pos = key
        contig = contig_of(chrom)
        if contig not in index.names:
            return [(v, None) for v in by_pos[key]]
        records = gf.with_retry(lambda: gr.search(url, index, chrom, pos, fetch, contig=contig))
        out = []
        for v in by_pos[key]:
            hit = None
            for r in records:
                hit = parse(r, v["ref"], v["alt"]) or hit
            out.append((v, hit))
        return out

    with ThreadPoolExecutor(threads) as pool:
        return [x for chunk in pool.map(one, sorted(by_pos)) for x in chunk]


def variants_to_score() -> list:
    """The rows annotated with gnomAD frequencies (study regions, new labels), SNVs only for AVI is
    decided by the caller."""
    with gzip.open(gf.VARIANTS, "rt", encoding="utf-8", newline="") as h:
        return [r for r in csv.DictReader(h, delimiter="\t") if gf.selected(r, "new")]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--tool", required=True, choices=("avi", "spliceai", "cadd"))
    ap.add_argument("--avi-file", type=Path)
    ap.add_argument("--threads", type=int, default=24)
    a = ap.parse_args(argv)
    variants = [v for v in variants_to_score() if len(v["ref"]) == 1 and len(v["alt"]) == 1]
    if a.tool == "avi":
        url, fetch = str(a.avi_file), local_fetch(a.avi_file)
        tbi = Path(str(a.avi_file) + ".tbi").read_bytes()
        with gzip.open(a.avi_file, "rt", encoding="utf-8") as h:
            header = next(line for line in h if not line.startswith("##")).rstrip("\n").split("\t")
        parse, columns = parse_avi(header), ("avi_raw", "avi_phred")
    else:
        url = SPLICEAI if a.tool == "spliceai" else CADD
        fetch = gr.fetch_range
        tbi = gf.with_retry(lambda: urllib.request.urlopen(url + ".tbi", timeout=120).read())
        parse = parse_spliceai if a.tool == "spliceai" else parse_cadd
        columns = ("spliceai_max",) if a.tool == "spliceai" else ("cadd_raw", "cadd_phred")
    index = gr.parse_tbi(tbi)
    chr_prefix = index.names[0].startswith("chr")
    contig_of = (lambda c: f"chr{c}") if chr_prefix else (lambda c: c)
    t0 = time.time()
    results = lookup_all(variants, url, index, contig_of, parse, fetch, a.threads)
    OUT.mkdir(parents=True, exist_ok=True)
    with gzip.open(OUT / f"{a.tool}.tsv.gz", "wt", encoding="utf-8", newline="") as h:
        w = csv.writer(h, delimiter="\t", lineterminator="\n")
        w.writerow(("variation_id", "chrom", "pos", "ref", "alt") + columns)
        for v, hit in sorted(results, key=lambda x: (x[0]["chrom"], int(x[0]["pos"]), x[0]["variation_id"])):
            w.writerow([v["variation_id"], v["chrom"], v["pos"], v["ref"], v["alt"]]
                       + [(hit or {}).get(c, "") for c in columns])
    source = {"tool": a.tool, "file": url, "index_sha256": hashlib.sha256(tbi).hexdigest(),
              "variants": len(results), "scored": sum(1 for _, h in results if h), "seconds": round(time.time() - t0)}
    (OUT / f"{a.tool}_source.json").write_text(json.dumps(source, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(source))
    return 0


if __name__ == "__main__":
    sys.exit(main())
