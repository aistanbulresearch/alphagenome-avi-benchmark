"""Conservation scores at the positions of the scored SNVs, from UCSC bigWig files by byte ranges.

    phyloP 100-way        comparator
    phastCons 470-way     AVI input
    Cactus 241-way phyloP AVI input (Zoonomia)

Each data block is fetched once for all positions it covers. Output:
data/derived/scores/conservation.tsv.gz (variation_id, chrom, pos, phyloP100, phastCons470,
phyloP241) and conservation_source.json.

Usage:
    python src/conservation_scores.py [--threads 16]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

import bigwig_remote as bw
import gnomad_frequencies as gf
import tool_scores as ts
from paths import DERIVED

UCSC = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/"
FILES = {"phyloP100": UCSC + "phyloP100way/hg38.phyloP100way.bw",
         "phastCons470": UCSC + "phastCons470way/hg38.phastCons470way.bw",
         "phyloP241": UCSC + "cactus241way/cactus241way.phyloP.bw"}
OUT = DERIVED / "scores" / "conservation.tsv.gz"


def values_for(url: str, positions: list, threads: int) -> dict:
    """{(chrom, pos): value or None} for 1-based positions on chromosomes named '1'..'22', 'X'."""
    f = bw.BigWig(url)
    by_block, out = defaultdict(list), {}
    for chrom, pos in positions:
        name = f"chr{chrom}"
        if name not in f.chroms:
            out[(chrom, pos)] = None
            continue
        blocks = gf.with_retry(lambda: f.blocks(f.chroms[name][0], pos - 1, pos))
        if not blocks:
            out[(chrom, pos)] = None
        for blk in blocks:
            by_block[blk].append((chrom, pos))

    def one(item):
        blk, keys = item
        vals = gf.with_retry(lambda: f.block_values(*blk))
        found = {}
        for chrom, pos in keys:
            for s, e, v in vals:
                if s <= pos - 1 < e:
                    found[(chrom, pos)] = float(v)
                    break
        return keys, found

    with ThreadPoolExecutor(threads) as pool:
        for keys, found in pool.map(one, by_block.items()):
            for k in keys:
                if k in found or k not in out:
                    out[k] = found.get(k)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--threads", type=int, default=16)
    a = ap.parse_args(argv)
    variants = [v for v in ts.variants_to_score() if len(v["ref"]) == 1 and len(v["alt"]) == 1]
    positions = sorted({(v["chrom"], int(v["pos"])) for v in variants})
    values, sources = {}, {}
    for name, url in FILES.items():
        values[name] = values_for(url, positions, a.threads)
        sources[name] = {"url": url, "positions": len(positions),
                         "with_value": sum(x is not None for x in values[name].values())}
        print(name, sources[name], flush=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(OUT, "wt", encoding="utf-8", newline="") as h:
        w = csv.writer(h, delimiter="\t", lineterminator="\n")
        w.writerow(("variation_id", "chrom", "pos") + tuple(FILES))
        for v in variants:
            key = (v["chrom"], int(v["pos"]))
            w.writerow([v["variation_id"], v["chrom"], v["pos"]]
                       + ["" if values[n].get(key) is None else f"{values[n][key]:g}" for n in FILES])
    (OUT.parent / "conservation_source.json").write_text(json.dumps(sources, indent=1) + "\n", encoding="utf-8",
                                                         newline="\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
