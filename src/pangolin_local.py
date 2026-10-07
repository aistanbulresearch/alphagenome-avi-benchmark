"""Pangolin run locally on the main set (new labels, untouched + validation chromosomes, rare, at least one
star, MANE-named) of the non-canonical intronic and UTR regions: every P/LP variant and a random B/LB
sample of 4,000 stratified by region (seed 20261007; list in data/derived/variants/pangolin_sample.tsv).

Follows pangolin.py of github.com/tkzeng/Pangolin (defaults: distance 50, masked): genes are GENCODE v38
`gene` features overlapping the variant; annotated splice sites are the exon ends of Ensembl_canonical
transcripts (the filter of the authors' create_db.py). Score per gene: the larger of the maximum gain
and the absolute maximum loss after masking (gains at annotated sites and losses at unannotated sites
set to 0); the variant's score is the maximum over its genes. The sequences of many variants are
predicted in one batch per model. A check against the package's compute_score on positions outside the
study sets runs first.

Runs in its own environment (Python 3.11, torch, pangolin without its VCF dependencies).

Output: data/derived/scores/pangolin.tsv.gz (variation_id, chrom, pos, ref, alt, pangolin_max, gene).

Usage:
    python src/pangolin_local.py --fasta <hg38.fa> --gtf <gencode.v38.annotation.gtf.gz> [--batch 16]
"""
from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import re
import time
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
VARIANTS = ROOT / "data" / "derived" / "variants" / "clinvar_variants.tsv.gz"
GNOMAD = ROOT / "data" / "derived" / "variants" / "gnomad_v4.1_joint.tsv.gz"
OUT = ROOT / "data" / "derived" / "scores" / "pangolin.tsv.gz"
SAMPLE = ROOT / "data" / "derived" / "variants" / "pangolin_sample.tsv"
PART2_TABLE = ROOT / "data" / "derived" / "part2" / "assay_variants.tsv.gz"
PART2_OUT = ROOT / "data" / "derived" / "scores" / "part2_pangolin.tsv.gz"
PART2_DATASETS = ("kitzman_pou1f1_ex2", "kitzman_wt1_ex9", "braun2018_mst1r_ex11", "cortes_lopez2022_cd19")
D = 50
FLANK = 5000
MAIN_REGIONS = ("donor_region", "acceptor_region", "near_intronic", "deep_intronic", "utr5", "utr3")
CODE = {"A": 1, "C": 2, "G": 3, "T": 4}
TABLE = np.asarray([[0, 0, 0, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], dtype=np.float32)


def one_hot(seq: str, strand: str) -> np.ndarray:
    """Pangolin's encoding: (4, L); the minus strand is reverse-complemented."""
    idx = np.frombuffer(bytes(CODE.get(c, 0) for c in seq.upper()), np.int8).astype(int)
    if strand == "-":
        idx = (5 - idx[::-1]) % 5
    return TABLE[idx].T


def read_annotation(gtf: Path) -> dict:
    """{chrom: (sorted gene starts, genes)}; gene = (start, end, strand, gene_id, sorted exon ends)."""
    genes, exons = {}, defaultdict(set)
    with gzip.open(gtf, "rt", encoding="utf-8") as h:
        for line in h:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] not in ("gene", "exon"):
                continue
            gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
            if f[2] == "gene":
                genes[gid] = (f[0], int(f[3]), int(f[4]), f[6])
            elif "Ensembl_canonical" in f[8]:
                exons[gid].update((int(f[3]), int(f[4])))
    by_chrom = defaultdict(list)
    for gid, (chrom, s, e, strand) in genes.items():
        by_chrom[chrom].append((s, e, strand, gid, np.array(sorted(exons.get(gid, ())))))
    return {c: ([g[0] for g in sorted(v)], sorted(v, key=lambda g: g[0])) for c, v in by_chrom.items()}


def genes_at(ann: dict, chrom: str, pos: int) -> list:
    starts, genes = ann.get(chrom, ([], []))
    i = bisect.bisect_right(starts, pos)
    return [g for g in genes[max(0, i - 3000):i] if g[0] <= pos <= g[1]]


def load_models():
    import torch
    from pkg_resources import resource_filename
    from pangolin.model import AR, L, W, Pangolin
    models = []
    for i in (0, 2, 4, 6):
        for j in (1, 2, 3):
            m = Pangolin(L, W, AR)
            m.load_state_dict(torch.load(resource_filename("pangolin", f"models/final.{j}.{i}.3.v2"),
                                         map_location="cpu", weights_only=True))
            m.eval()
            models.append(m)
    return models


def predict(models: list, x: np.ndarray) -> np.ndarray:
    """(n, 4, 2D+1): mean over the three models of each tissue group of the splice-usage track."""
    import torch
    t = torch.from_numpy(x)
    out = np.empty((x.shape[0], 4, 2 * D + 1), dtype=np.float32)
    with torch.no_grad():
        for j in range(4):
            out[:, j] = np.mean([m(t)[:, [1, 4, 7, 10][j], :].numpy() for m in models[3 * j:3 * j + 3]], axis=0)
    return out


def masked_score(diff: np.ndarray, sites: np.ndarray, pos: int) -> float:
    """diff: (4, 2D+1) alt - ref in genomic order. Pangolin's loss/gain, masking, maximum magnitude."""
    loss = diff[np.argmin(diff, axis=0), np.arange(diff.shape[1])].copy()
    gain = diff[np.argmax(diff, axis=0), np.arange(diff.shape[1])].copy()
    p = sites - (pos - D)
    p = p[(p >= 0) & (p < len(loss))]
    if len(p):
        gain[p] = np.minimum(gain[p], 0)
        other = ~np.isin(np.arange(len(loss)), p)
        loss[other] = np.maximum(loss[other], 0)
    else:
        loss[:] = np.maximum(loss, 0)
    return float(max(gain.max(), -loss.min()))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fasta", required=True, type=Path)
    ap.add_argument("--gtf", required=True, type=Path)
    ap.add_argument("--batch", type=int, default=16)
    ap.add_argument("--part2", action="store_true", help="the SNVs of the four single-exon Part 2 datasets")
    a = ap.parse_args(argv)
    from pyfaidx import Fasta
    import pangolin_check
    fa, ann, models = Fasta(str(a.fasta)), read_annotation(a.gtf), load_models()
    if a.part2:
        main_rows = rows = pangolin_check.part2_rows(PART2_TABLE, PART2_DATASETS)
        out, id_column = PART2_OUT, "allele"
        print("part 2 alleles", len(rows), flush=True)
    else:
        main_rows = pangolin_check.main_set_rows(VARIANTS, GNOMAD, MAIN_REGIONS)
        rows = pangolin_check.sample(main_rows)
        out, id_column = OUT, "variation_id"
        with open(SAMPLE, "w", encoding="utf-8", newline="") as h:
            w = csv.writer(h, delimiter="\t", lineterminator="\n")
            w.writerow(("variation_id", "label", "region"))
            w.writerows((r["variation_id"], r["label"], r["region"]) for r in rows)
        print("main set", len(main_rows), "sample", len(rows), flush=True)
    pangolin_check.self_check(fa, ann, models, main_rows)
    print("self-check passed", flush=True)

    best = {}

    def run(part: list) -> None:
        x = np.stack([j[2] for j in part] + [j[3] for j in part])
        y = predict(models, x)
        n = len(part)
        for i, (vid, strand, _, _, genes, pos) in enumerate(part):
            ref_y, alt_y = y[i], y[n + i]
            if strand == "-":
                ref_y, alt_y = ref_y[:, ::-1], alt_y[:, ::-1]
            for g in genes:
                s = masked_score(alt_y - ref_y, g[4], pos)
                if vid not in best or s > best[vid][0]:
                    best[vid] = (s, g[3])

    pending, started = [], time.time()       # (variation_id, strand, ref one-hot, alt one-hot, genes, pos)
    for n_row, r in enumerate(rows, 1):
        chrom, pos, ref, alt = "chr" + r["chrom"], int(r["pos"]), r["ref"], r["alt"]
        seq = fa[chrom][pos - FLANK - D - 1:pos + FLANK + D].seq
        if len(seq) == 2 * (FLANK + D) + 1 and seq[FLANK + D].upper() == ref:
            alt_seq = seq[:FLANK + D] + alt + seq[FLANK + D + 1:]
            genes = genes_at(ann, chrom, pos)
            for strand in sorted({g[2] for g in genes}):
                sites = [g for g in genes if g[2] == strand]
                pending.append((r["variation_id"], strand, one_hot(seq, strand), one_hot(alt_seq, strand), sites, pos))
        if len(pending) >= a.batch or (n_row == len(rows) and pending):
            run(pending)
            pending = []
        if n_row % 500 == 0 or n_row == len(rows):
            print(f"{n_row}/{len(rows)} {time.time() - started:.0f} s", flush=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out, "wt", encoding="utf-8", newline="") as h:
        w = csv.writer(h, delimiter="\t", lineterminator="\n")
        w.writerow((id_column, "chrom", "pos", "ref", "alt", "pangolin_max", "gene"))
        for r in rows:
            s = best.get(r["variation_id"])
            w.writerow((r["variation_id"], r["chrom"], r["pos"], r["ref"], r["alt"],
                        "" if s is None else f"{s[0]:.4f}", "" if s is None else s[1]))
    print("variants", len(rows), "scored", len(best))
    return 0


if __name__ == "__main__":
    sys.exit(main())
