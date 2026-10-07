"""Scores for the Part 2 assay SNVs: AVI (Atlas API), SpliceAI raw and masked, CADD and conservation, with
the readers of Part 1. Every distinct allele (chrom, pos, ref, alt) of the main and comparison rows is
scored once; the analysis joins scores back by allele.

Output: data/derived/scores/part2_<tool>.tsv.gz (allele, chrom, pos, ref, alt, score columns) and
part2_<tool>_source.json.

Usage:
    python src/part2_scores.py --tool {avi,spliceai,spliceai_masked,cadd,conservation} [--threads 16]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import sys
import time
import urllib.request

import pandas as pd

import gnomad_frequencies as gf
import gnomad_remote as gr
import tool_scores as ts
from paths import DERIVED

TABLE = DERIVED / "part2" / "assay_variants.tsv.gz"
OUT = DERIVED / "scores"


def alleles() -> list:
    """Distinct SNV alleles of the main and comparison rows, as the dicts the Part 1 readers expect."""
    d = pd.read_csv(TABLE, sep="\t", dtype=str, keep_default_na=False)
    d = d[d["use"].isin(["main", "comparison"]) & (d["variant_type"] == "SNV")]
    d = d.drop_duplicates(["chrom", "pos", "ref", "alt"])
    return [{"variation_id": f"{r.chrom}:{r.pos}:{r.ref}:{r.alt}", "chrom": r.chrom, "pos": r.pos,
             "ref": r.ref, "alt": r.alt} for r in d.itertuples()]


def write(tool: str, columns: tuple, results: list, source: dict) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with gzip.open(OUT / f"part2_{tool}.tsv.gz", "wt", encoding="utf-8", newline="") as h:
        w = csv.writer(h, delimiter="\t", lineterminator="\n")
        w.writerow(("allele", "chrom", "pos", "ref", "alt") + columns)
        for v, hit in results:
            w.writerow([v["variation_id"], v["chrom"], v["pos"], v["ref"], v["alt"]]
                       + [(hit or {}).get(c, "") for c in columns])
    source = {**source, "variants": len(results), "scored": sum(1 for _, h in results if h)}
    (OUT / f"part2_{tool}_source.json").write_text(json.dumps(source, indent=1) + "\n", encoding="utf-8",
                                                   newline="\n")
    print(json.dumps(source), flush=True)


def tabix_tool(tool: str, variants: list, threads: int) -> None:
    url = ts.SPLICEAI.get(tool, ts.CADD)
    tbi = gf.with_retry(lambda: urllib.request.urlopen(url + ".tbi", timeout=120).read())
    index = gr.parse_tbi(tbi)
    contig_of = (lambda c: f"chr{c}") if index.names[0].startswith("chr") else (lambda c: c)
    parse = ts.parse_spliceai if tool in ts.SPLICEAI else ts.parse_cadd
    columns = ("spliceai_max",) if tool in ts.SPLICEAI else ("cadd_raw", "cadd_phred")
    t0 = time.time()
    results = ts.lookup_all(variants, url, index, contig_of, parse, gr.fetch_range, threads)
    write(tool, columns, results, {"tool": tool, "file": url, "seconds": round(time.time() - t0)})


def avi(variants: list, workers: int) -> None:
    import avi_atlas as aa
    from alphagenome.atlas import atlas
    from alphagenome.data import genome
    client = atlas.create(os.environ["ALPHAGENOME_API_KEY"])

    def query(ids):
        return client.query_variants([genome.Variant.from_str(i) for i in ids], requested_scorers=[aa.SCORER],
                                     progress_bar=False, max_workers=workers)
    results = []
    for k in range(0, len(variants), 500):
        chunk = variants[k:k + 500]
        for v, raw, q, status in aa.run_chunk(query, chunk):
            results.append((v, {"avi_raw": raw, "avi_quantile": q, "status": status} if status == "ok" else None))
        if (k // 500) % 20 == 0:
            print(f"{min(k + 500, len(variants))}/{len(variants)}", flush=True)
    write("avi", ("avi_raw", "avi_quantile"), results, {"tool": "AVI", "scorer": aa.SCORER})


def conservation(variants: list, threads: int) -> None:
    import conservation_scores as cs
    positions = sorted({(v["chrom"], int(v["pos"])) for v in variants})
    values = {n: cs.values_for(url, positions, threads) for n, url in cs.FILES.items()}
    results = []
    for v in variants:
        key = (v["chrom"], int(v["pos"]))
        hit = {n: f"{values[n][key]:g}" for n in cs.FILES if values[n].get(key) is not None}
        results.append((v, hit or None))
    write("conservation", tuple(cs.FILES), results, {"tool": "conservation", "files": cs.FILES})


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--tool", required=True, choices=("avi", "spliceai", "spliceai_masked", "cadd", "conservation"))
    ap.add_argument("--threads", type=int, default=16)
    a = ap.parse_args(argv)
    variants = alleles()
    print("alleles", len(variants), flush=True)
    if a.tool == "avi":
        avi(variants, 20)
    elif a.tool == "conservation":
        conservation(variants, a.threads)
    else:
        tabix_tool(a.tool, variants, a.threads)
    return 0


if __name__ == "__main__":
    sys.exit(main())
