"""gnomAD v4.1 joint frequencies for the variants of the ClinVar variant list.

Reads `data/derived/variants/clinvar_variants.tsv.gz`, selects the rows to annotate, looks each position
up in the gnomAD v4.1 joint sites VCF of its chromosome through the remote tabix reader
(src/gnomad_remote.py) and writes:

    data/derived/variants/gnomad_v4.1_joint.tsv.gz   one row per annotated variant
    data/derived/variants/gnomad_v4.1_joint_sources.json   per chromosome: URL, index SHA-256, ETag

Each chromosome is first written to its own file under data/derived/variants/gnomad_parts/ (with its
source record); a rerun skips finished chromosomes, so an interrupted run resumes.

Alleles are compared after the same trimming that the ClinVar keys use (src/clinvar_records.py).

Usage:
    python src/gnomad_frequencies.py [--rows new|all] [--threads 24]
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

import clinvar_records as cr
import gnomad_remote as gr
from paths import DERIVED

VARIANTS = DERIVED / "variants" / "clinvar_variants.tsv.gz"
OUT = DERIVED / "variants" / "gnomad_v4.1_joint.tsv.gz"
SOURCES = DERIVED / "variants" / "gnomad_v4.1_joint_sources.json"
PARTS = DERIVED / "variants" / "gnomad_parts"
CHROMS = [str(i) for i in range(1, 23)] + ["X"]
COLUMNS = ("variation_id", "chrom", "pos", "ref", "alt", "in_gnomad", "filter") + gr.FIELDS


def selected(row: dict, rows: str) -> bool:
    """Rows to annotate: chromosomes 1-22 and X, at least one star, the regions of the study (intronic,
    UTR, synonymous); 'new' keeps the new-label set."""
    if row["chrom"] not in CHROMS or not row["stars"].isdigit() or int(row["stars"]) < 1:
        return False
    if row["family"] not in ("intronic", "utr") and row["region"] != "synonymous":
        return False
    return rows == "all" or row["new_label"] == "yes"


def with_retry(call, attempts: int = 16):
    """Retry network errors and truncated responses with back-off (up to about 12 minutes in all)."""
    for i in range(attempts):
        try:
            return call()
        except (OSError, ValueError):
            if i == attempts - 1:
                raise
            time.sleep(min(60, 2 ** i))


def annotate_chromosome(chrom: str, variants: list, threads: int) -> tuple[list, dict]:
    url = gr.vcf_url(chrom)

    def get_index():
        with urllib.request.urlopen(url + ".tbi", timeout=120) as resp:
            return resp.read()

    def get_head():
        with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=120) as resp:
            return resp.headers.get("ETag"), resp.headers.get("Last-Modified")

    tbi = with_retry(get_index)
    etag, modified = with_retry(get_head)
    index = gr.parse_tbi(tbi)
    by_pos = defaultdict(list)
    for v in variants:
        by_pos[int(v["pos"])].append(v)

    def one(pos: int) -> list:
        records = with_retry(lambda: gr.search(url, index, chrom, pos))
        found = {}
        for r in records:
            key = cr.normalize_allele(chrom, r[1], r[3], r[4])
            if key is not None:
                info = gr.parse_info(r[7])
                info["FILTER"] = r[6]
                found[key] = info
        out = []
        for v in by_pos[pos]:
            info = found.get(f"{chrom}:{v['pos']}:{v['ref']}:{v['alt']}")
            out.append({"variation_id": v["variation_id"], "chrom": chrom, "pos": v["pos"], "ref": v["ref"],
                        "alt": v["alt"], "in_gnomad": "yes" if info else "no",
                        "filter": info["FILTER"] if info else "",
                        **{f: (info or {}).get(f, "") for f in gr.FIELDS}})
        return out

    with ThreadPoolExecutor(threads) as pool:
        rows = [r for chunk in pool.map(one, sorted(by_pos)) for r in chunk]
    source = {"url": url, "index_sha256": hashlib.sha256(tbi).hexdigest(), "etag": etag,
              "last_modified": modified, "positions": len(by_pos), "variants": len(rows)}
    return rows, source


def merge(chroms: list, rows: str) -> None:
    """Concatenate the per-chromosome files into OUT and their source records into SOURCES."""
    sources = {}
    with gzip.open(OUT, "wt", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for chrom in chroms:
            with gzip.open(PARTS / f"chr{chrom}.tsv.gz", "rt", encoding="utf-8", newline="") as h:
                writer.writerows(csv.DictReader(h, delimiter="\t"))
            sources[chrom] = json.loads((PARTS / f"chr{chrom}.json").read_text(encoding="utf-8"))
    SOURCES.write_text(json.dumps({"rows": rows, "chromosomes": sources}, indent=1) + "\n",
                       encoding="utf-8", newline="\n")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--rows", choices=("new", "all"), default="new")
    ap.add_argument("--threads", type=int, default=24)
    a = ap.parse_args(argv)
    per_chrom = defaultdict(list)
    with gzip.open(VARIANTS, "rt", encoding="utf-8", newline="") as h:
        for row in csv.DictReader(h, delimiter="\t"):
            if selected(row, a.rows):
                per_chrom[row["chrom"]].append(row)
    PARTS.mkdir(parents=True, exist_ok=True)
    for chrom in CHROMS:
        part, meta = PARTS / f"chr{chrom}.tsv.gz", PARTS / f"chr{chrom}.json"
        if not per_chrom[chrom] or (part.exists() and meta.exists()):
            continue
        t0 = time.time()
        rows, source = annotate_chromosome(chrom, per_chrom[chrom], a.threads)
        tmp = part.with_suffix(".tmp")
        with gzip.open(tmp, "wt", encoding="utf-8", newline="") as h:
            writer = csv.DictWriter(h, COLUMNS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(sorted(rows, key=lambda r: (int(r["pos"]), r["variation_id"])))
        tmp.replace(part)
        meta.write_text(json.dumps(source, indent=1) + "\n", encoding="utf-8", newline="\n")
        print(f"chr{chrom}\t{len(rows)} variants\t{time.time() - t0:.0f} s", flush=True)
    merge([c for c in CHROMS if per_chrom[c]], a.rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
