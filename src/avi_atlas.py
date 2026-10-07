"""AVI scores for the SNVs of the ClinVar variant list from the AlphaGenome Atlas API.

The API key is read from the environment variable ALPHAGENOME_API_KEY and is never written anywhere.
Scorer `AVI_SCORE`: the raw score and its 'quantiles' layer (the calibrated score the API returns).
Variants are sent in chunks; a chunk that fails is retried and then sent one variant at a time, so a
variant the API rejects (for example a reference base that does not match) is recorded, not fatal.
Finished chunks are appended to a partial file, so an interrupted run resumes.

Output: data/derived/scores/avi.tsv.gz and avi_source.json.

Usage:
    python src/avi_atlas.py [--chunk 500] [--workers 20]
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import json
import os
import sys
import time
from pathlib import Path
from typing import Callable

import tool_scores as ts
from paths import DERIVED

OUT = DERIVED / "scores"
SCORER = "AVI_SCORE"
COLUMNS = ("variation_id", "chrom", "pos", "ref", "alt", "avi_raw", "avi_quantile", "status")


def variant_id(v: dict) -> str:
    return f"chr{v['chrom']}:{v['pos']}:{v['ref']}>{v['alt']}"


def scores_of(result) -> dict:
    """{variant id: (raw, quantile)} from the AnnData the client returns for AVI_SCORE."""
    ad = result[SCORER]
    quant = ad.layers.get("quantiles") if hasattr(ad.layers, "get") else None
    out = {}
    for i, vid in enumerate(ad.obs["variant"].astype(str)):
        raw = float(ad.X[i, 0])
        q = float(quant[i, 0]) if quant is not None else float("nan")
        out[vid] = (raw, q)
    return out


def run_chunk(query: Callable, chunk: list, retries: int = 4) -> list:
    """[(variant, raw, quantile, status)] for one chunk; per-variant fallback after repeated failure."""
    ids = [variant_id(v) for v in chunk]
    for attempt in range(retries):
        try:
            got = scores_of(query(ids))
            return [(v, *got.get(i, ("", "")), "ok" if i in got else "not returned") for v, i in zip(chunk, ids)]
        except Exception:                                        # noqa: BLE001  network or a rejected variant
            time.sleep(min(60, 2 ** attempt))
    out = []
    for v, i in zip(chunk, ids):
        try:
            got = scores_of(query([i]))
            out.append((v, *got.get(i, ("", "")), "ok" if i in got else "not returned"))
        except Exception as err:                                 # noqa: BLE001
            out.append((v, "", "", f"error: {type(err).__name__}"))
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--chunk", type=int, default=500)
    ap.add_argument("--workers", type=int, default=20)
    a = ap.parse_args(argv)
    from alphagenome.atlas import atlas                         # imported here: optional dependency
    from alphagenome.data import genome
    client = atlas.create(os.environ["ALPHAGENOME_API_KEY"])

    def query(ids):
        return client.query_variants([genome.Variant.from_str(i) for i in ids], requested_scorers=[SCORER],
                                     progress_bar=False, max_workers=a.workers)

    variants = [v for v in ts.variants_to_score() if len(v["ref"]) == 1 and len(v["alt"]) == 1]
    OUT.mkdir(parents=True, exist_ok=True)
    partial = OUT / "avi.partial.tsv"
    done = set()
    if partial.exists():
        with open(partial, encoding="utf-8", newline="") as h:
            done = {r["variation_id"] for r in csv.DictReader(h, delimiter="\t")}
    todo = [v for v in variants if v["variation_id"] not in done]
    started = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    with open(partial, "a", encoding="utf-8", newline="") as h:
        w = csv.writer(h, delimiter="\t", lineterminator="\n")
        if not done:
            w.writerow(COLUMNS)
        for k in range(0, len(todo), a.chunk):
            for v, raw, q, status in run_chunk(query, todo[k:k + a.chunk]):
                w.writerow([v["variation_id"], v["chrom"], v["pos"], v["ref"], v["alt"], raw, q, status])
            h.flush()
            print(f"{min(k + a.chunk, len(todo))}/{len(todo)}", flush=True)
    with open(partial, encoding="utf-8", newline="") as h:
        rows = list(csv.DictReader(h, delimiter="\t"))
    with gzip.open(OUT / "avi.tsv.gz", "wt", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["chrom"], int(r["pos"]), r["variation_id"])))
    source = {"tool": "AVI", "api": "AlphaGenome Atlas (alphagenome client)", "scorer": SCORER,
              "client_version": _version(), "started_utc": started,
              "finished_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
              "variants": len(rows), "scored": sum(r["status"] == "ok" for r in rows)}
    (OUT / "avi_source.json").write_text(json.dumps(source, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(source))
    return 0


def _version() -> str:
    try:
        from importlib.metadata import version
        return version("alphagenome")
    except Exception:                                            # noqa: BLE001
        return "unknown"


if __name__ == "__main__":
    sys.exit(main())
