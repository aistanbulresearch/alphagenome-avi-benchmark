"""Helpers for src/pangolin_local.py: the rows of the main set, and a check of the batched scoring against
the package's own compute_score on positions outside the study sets."""
from __future__ import annotations

import csv
import gzip
import sys
import types
from pathlib import Path

import numpy as np

RARE = 0.001


def main_set_rows(variants: Path, gnomad: Path, regions: tuple) -> list:
    """New labels, SNVs, >=1 star, MANE-named, untouched or validation chromosomes, rare in gnomAD v4.1."""
    freq = {}
    with gzip.open(gnomad, "rt", encoding="utf-8", newline="") as h:
        for g in csv.DictReader(h, delimiter="\t"):
            x = g["fafmax_faf95_max_joint"] or g["AF_joint"]
            freq[g["variation_id"]] = (g["in_gnomad"] == "yes", float(x) if x not in ("", ".") else None)
    out = []
    with gzip.open(variants, "rt", encoding="utf-8", newline="") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if not (r["new_label"] == "yes" and r["variant_type"] == "snv" and r["region"] in regions
                    and r["stars"].isdigit() and int(r["stars"]) >= 1 and r["mane_select"] == "yes"
                    and r["chromosome_group"] in ("untouched", "validation")):
                continue
            present, f = freq.get(r["variation_id"], (False, None))
            if not present or f is None or f < RARE:
                out.append(r)
    return out


def self_check(fa, ann, models, rows: list) -> None:
    """Batched loss/gain must equal pangolin.pangolin.compute_score on 4 positions outside the study."""
    import pangolin_local as pl
    for name in ("vcf", "gffutils", "pyfastx"):                 # unused by compute_score
        sys.modules.setdefault(name, types.ModuleType(name))
    from pangolin.pangolin import compute_score
    study = {("chr" + r["chrom"], int(r["pos"])) for r in rows}
    checks = [("chr17", 43_100_000 + 41 * k) for k in range(12)]
    checks = [c for c in checks if c not in study][:4]
    for chrom, pos in checks:
        seq = fa[chrom][pos - pl.FLANK - pl.D - 1:pos + pl.FLANK + pl.D].seq
        ref = seq[pl.FLANK + pl.D].upper()
        alt = "A" if ref != "A" else "G"
        alt_seq = seq[:pl.FLANK + pl.D] + alt + seq[pl.FLANK + pl.D + 1:]
        for strand in ("+", "-"):
            loss_p, gain_p = compute_score(seq, alt_seq, strand, pl.D, models)
            y = pl.predict(models, np.stack([pl.one_hot(seq, strand), pl.one_hot(alt_seq, strand)]))
            r_, a_ = (y[0], y[1]) if strand == "+" else (y[0][:, ::-1], y[1][:, ::-1])
            diff = a_ - r_
            loss = diff[np.argmin(diff, axis=0), np.arange(diff.shape[1])]
            gain = diff[np.argmax(diff, axis=0), np.arange(diff.shape[1])]
            if np.abs(loss - loss_p).max() > 1e-4 or np.abs(gain - gain_p).max() > 1e-4:
                raise SystemExit(f"self-check failed at {chrom}:{pos} strand {strand}")
