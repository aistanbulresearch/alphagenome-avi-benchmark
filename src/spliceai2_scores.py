"""SpliceAI2 scores (precomputed v2.0, GRCh38) for every Part 1 SNV row and every Part 2 SNV allele, read
by byte ranges from the per-chromosome files (src/unindexed_remote.py); STUDY_DESIGN.md section 4.4.

Score: `spliceai2_summary_score`, the maximum over the gene/strand rows of the allele. The largest of the
four strongest delta scores (`*_delta_score_0`) is kept beside it, and the number of rows where the two
differ is recorded. Access needs the user's Hugging Face token (HF_TOKEN) with the dataset terms accepted;
the scores are for non-commercial research use and are not redistributed (data/derived is not in git).

Output: data/derived/scores/spliceai2.tsv.gz (Part 1, by variation_id), part2_spliceai2.tsv.gz (Part 2,
by allele), spliceai2_source.json; per-chromosome parts in spliceai2_parts/ make the run resumable.

Usage:
    python src/spliceai2_scores.py [--threads 32] [--step 8000000] [--chroms 21 ...]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import sys
import threading
import time
import urllib.error
import urllib.request

import gnomad_frequencies as gf
import part2_scores as p2
import tool_scores as ts
import unindexed_remote as ur
from paths import DERIVED

REPO = "https://huggingface.co/datasets/illumina-ai/SpliceAI2-data"
PREFIX = "precomputed_scores_v2.0/GRCh38/snv/"
OUT = DERIVED / "scores"
PARTS = OUT / "spliceai2_parts"
DELTA0 = ("donor_gain_delta_score_0", "donor_loss_delta_score_0", "acceptor_gain_delta_score_0",
          "acceptor_loss_delta_score_0")
FIELDS = ("spliceai2_summary", "spliceai2_maxdelta", "rows", "differs")


def columns(header: list) -> dict:
    return {"summary": header.index("spliceai2_summary_score"), "delta0": [header.index(c) for c in DELTA0],
            "width": len(header)}


def combine(rows: list, ref: str, alt: str, cols: dict) -> dict | None:
    """One score per allele from its gene/strand rows, or None if the allele has no row."""
    best, best_delta, n, differs = None, None, 0, 0
    for r in rows:
        if r[2] != ref or r[3] != alt:
            continue
        summary = float(r[cols["summary"]])
        delta = max(float(r[i]) for i in cols["delta0"])
        n += 1
        differs += abs(summary - delta) > 1e-6
        best = summary if best is None else max(best, summary)
        best_delta = delta if best_delta is None else max(best_delta, delta)
    if n == 0:
        return None
    return {"spliceai2_summary": f"{best:g}", "spliceai2_maxdelta": f"{best_delta:g}", "rows": n,
            "differs": differs}


class Remote:
    """Range reads through the signed download address of each file, renewed every 20 minutes or on 403."""

    def __init__(self, token: str):
        self.token, self.signed, self.lock = token, {}, threading.Lock()

    def resolve(self, url: str) -> str:
        req = urllib.request.Request(url, headers={"Authorization": f"Bearer {self.token}", "Range": "bytes=0-0"})
        with urllib.request.urlopen(req, timeout=120) as r:
            return r.geturl()

    def address(self, url: str, renew: bool = False) -> str:
        with self.lock:
            hit = self.signed.get(url)
            if renew or hit is None or time.time() - hit[1] > 1200:
                self.signed[url] = (self.resolve(url), time.time())
            return self.signed[url][0]

    def fetch(self, url: str, start: int, end: int) -> bytes:
        def once():
            try:
                req = urllib.request.Request(self.address(url), headers={"Range": f"bytes={start}-{end - 1}"})
                with urllib.request.urlopen(req, timeout=120) as r:
                    return r.read()
            except urllib.error.HTTPError as e:
                if e.code in (401, 403):
                    self.address(url, renew=True)
                raise
        return gf.with_retry(once)


def file_sizes(token: str) -> dict:
    req = urllib.request.Request(f"https://huggingface.co/api/datasets/illumina-ai/SpliceAI2-data/tree/main/"
                                 f"{PREFIX.rstrip('/')}", headers={"Authorization": f"Bearer {token}"})
    with urllib.request.urlopen(req, timeout=120) as r:
        listing = json.load(r)
    return {x["path"].split("/")[-1]: x.get("size") for x in listing if x.get("type") == "file"}


def header_of(remote: Remote, url: str) -> list:
    data = remote.fetch(url, 0, ur.PROBE)
    text = b"".join(b for _, b in ur.gr.inflate_blocks(data))
    return text.split(b"\n", 1)[0].decode().split("\t")


def run_chromosome(remote: Remote, chrom: str, size: int, keys: set, cols: dict, step: int, threads: int) -> None:
    url = f"{REPO}/resolve/main/{PREFIX}chr{chrom}.tsv.gz"
    name = f"chr{chrom}"
    t0 = time.time()
    smap = ur.sparse_map(url, size, name, remote.fetch, step, threads=threads)
    found = ur.fetch_positions(url, size, smap, name, [p for p, _, _ in keys], remote.fetch, threads=threads)
    PARTS.mkdir(parents=True, exist_ok=True)
    tmp = PARTS / f"chr{chrom}.tsv.gz.tmp"
    with gzip.open(tmp, "wt", encoding="utf-8", newline="") as h:
        w = csv.writer(h, delimiter="\t", lineterminator="\n")
        w.writerow(("pos", "ref", "alt") + FIELDS)
        for pos, ref, alt in sorted(keys):
            hit = combine(found.get(pos, []), ref, alt, cols)
            w.writerow([pos, ref, alt] + [(hit or {}).get(c, "") for c in FIELDS])
    tmp.replace(PARTS / f"chr{chrom}.tsv.gz")
    print(f"chr{chrom}: {len(keys)} alleles, {len(found)} positions found, map {len(smap)}, "
          f"{time.time() - t0:.0f} s", flush=True)


def read_parts() -> dict:
    out = {}
    for path in sorted(PARTS.glob("chr*.tsv.gz")):
        chrom = path.name[3:].split(".")[0]
        with gzip.open(path, "rt", encoding="utf-8", newline="") as h:
            for r in csv.DictReader(h, delimiter="\t"):
                out[(chrom, int(r["pos"]), r["ref"], r["alt"])] = r
    return out


def write(path, key_name: str, variants: list, scores: dict) -> int:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as h:
        w = csv.writer(h, delimiter="\t", lineterminator="\n")
        w.writerow((key_name, "chrom", "pos", "ref", "alt") + FIELDS)
        scored = 0
        for v in variants:
            r = scores.get((v["chrom"], int(v["pos"]), v["ref"], v["alt"]), {})
            scored += bool(r.get("spliceai2_summary"))
            w.writerow([v["variation_id"], v["chrom"], v["pos"], v["ref"], v["alt"]] + [r.get(c, "") for c in FIELDS])
    return scored


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--step", type=int, default=8_000_000)
    ap.add_argument("--chroms", nargs="*", help="only these chromosomes (a trial run); parts are kept")
    a = ap.parse_args(argv)
    token = os.environ["HF_TOKEN"]
    remote = Remote(token)
    part1 = [v for v in ts.variants_to_score() if len(v["ref"]) == 1 and len(v["alt"]) == 1]
    part2 = p2.alleles()
    keys: dict = {}
    for v in part1 + part2:
        keys.setdefault(v["chrom"], set()).add((int(v["pos"]), v["ref"], v["alt"]))
    sizes = file_sizes(token)
    first_url = f"{REPO}/resolve/main/{PREFIX}chr21.tsv.gz"
    header = header_of(remote, first_url)
    cols = columns(header)
    print(f"part 1 rows {len(part1)}, part 2 alleles {len(part2)}, chromosomes {len(keys)}", flush=True)
    for chrom in sorted(keys, key=lambda c: -sizes.get(f"chr{c}.tsv.gz", 0) or 0):
        if (PARTS / f"chr{chrom}.tsv.gz").exists() or (a.chroms and chrom not in a.chroms):
            continue
        size = sizes.get(f"chr{chrom}.tsv.gz")
        if not size:
            print(f"chr{chrom}: no file", flush=True)
            continue
        if header_of(remote, f"{REPO}/resolve/main/{PREFIX}chr{chrom}.tsv.gz") != header:
            raise ValueError(f"chr{chrom}: header differs from chr21")
        run_chromosome(remote, chrom, size, keys[chrom], cols, a.step, a.threads)
    if a.chroms:
        return 0
    scores = read_parts()
    n1 = write(OUT / "spliceai2.tsv.gz", "variation_id", part1, scores)
    n2 = write(OUT / "part2_spliceai2.tsv.gz", "allele", part2, scores)
    source = {"tool": "spliceai2", "files": f"{REPO} {PREFIX}chr*.tsv.gz", "version": "precomputed_scores_v2.0",
              "header_sha256": hashlib.sha256("\t".join(header).encode()).hexdigest(),
              "part1_rows": len(part1), "part1_scored": n1, "part2_alleles": len(part2), "part2_scored": n2,
              "rows_where_summary_differs_from_max_delta0": sum(int(r["differs"] or 0) for r in scores.values()),
              "rows_read": sum(int(r["rows"] or 0) for r in scores.values())}
    (OUT / "spliceai2_source.json").write_text(json.dumps(source, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(source), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
