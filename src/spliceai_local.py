"""SpliceAI 1.3.1 run locally with a +/-4,999-base window for the deep intronic SNVs of the new-label set.

This is the score form of the ClinGen thresholds (raw maximum delta score, distance 4,999; Walker et
al. 2023). The precomputed files use distance 50, which can miss pseudoexon-creating variants whose
partner splice site lies further away.

Runs in its own environment (Python 3.11 with spliceai 1.3.1, tensorflow and pyfaidx; pysam is not
needed). The scoring follows spliceai.utils.get_delta_scores with mask = 0 for single-nucleotide
variants, with the sequences of many variants predicted in one batch per model. A check against the
package function on positions outside the study sets runs first and stops the run on any difference.

Output: data/derived/scores/spliceai_d4999.tsv.gz (variation_id, chrom, pos, ref, alt, spliceai_max,
DS_AG, DS_AL, DS_DG, DS_DL, gene).

Usage (in the SpliceAI environment):
    python src/spliceai_local.py --fasta <hg38.fa> [--batch 32]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
VARIANTS = ROOT / "data" / "derived" / "variants" / "clinvar_variants.tsv.gz"
OUT = ROOT / "data" / "derived" / "scores" / "spliceai_d4999.tsv.gz"
DIST = 4999
CODES = {"A": 1, "C": 2, "G": 3, "T": 4}
TABLE = np.asarray([[0, 0, 0, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], dtype=np.float32)


def one_hot(seq: str) -> np.ndarray:
    """spliceai's encoding (A, C, G, T; anything else is all zeros), without np.fromstring."""
    return TABLE[np.frombuffer(bytes(CODES.get(c, 0) for c in seq.upper()), np.int8)]


def inputs(record, ann, dist: int = DIST) -> list:
    """[(gene, strand, x_ref, x_alt)] for one SNV, as spliceai.utils.get_delta_scores builds them."""
    cov = 2 * dist + 1
    wid = 10000 + cov
    genes, strands, idxs = ann.get_name_and_strand(record.chrom, record.pos)
    if len(idxs) == 0:
        return []
    chrom = ("chr" + record.chrom) if not record.chrom.startswith("chr") else record.chrom
    seq = ann.ref_fasta[chrom][record.pos - wid // 2 - 1:record.pos + wid // 2].seq
    if len(seq) != wid or seq[wid // 2:wid // 2 + 1].upper() != record.ref:
        return []
    out = []
    for i in range(len(idxs)):
        d0, d1, _ = ann.get_pos_data(idxs[i], record.pos)
        pad = [max(wid // 2 + d0, 0), max(wid // 2 - d1, 0)]
        x_ref = "N" * pad[0] + seq[pad[0]:wid - pad[1]] + "N" * pad[1]
        x_alt = x_ref[:wid // 2] + record.alt + x_ref[wid // 2 + 1:]
        r, a = one_hot(x_ref), one_hot(x_alt)
        if strands[i] == "-":
            r, a = r[::-1, ::-1], a[::-1, ::-1]
        out.append((genes[i], strands[i], r, a))
    return out


def deltas(y_ref: np.ndarray, y_alt: np.ndarray, strand: str) -> tuple:
    """(DS_AG, DS_AL, DS_DG, DS_DL) for one gene, unmasked."""
    if strand == "-":
        y_ref, y_alt = y_ref[::-1], y_alt[::-1]
    ag = (y_alt[:, 1] - y_ref[:, 1]).max()
    al = (y_ref[:, 1] - y_alt[:, 1]).max()
    dg = (y_alt[:, 2] - y_ref[:, 2]).max()
    dl = (y_ref[:, 2] - y_alt[:, 2]).max()
    return float(ag), float(al), float(dg), float(dl)


def score_batch(items: list, models: list) -> list:
    """items: [(key, gene, strand, x_ref, x_alt)] -> [(key, gene, (ag, al, dg, dl))]."""
    x = np.stack([it[3] for it in items] + [it[4] for it in items])
    y = np.mean([m.predict(x, verbose=0) for m in models], axis=0)
    n = len(items)
    return [(it[0], it[1], deltas(y[k], y[n + k], it[2])) for k, it in enumerate(items)]


def best(per_gene: list) -> tuple:
    """The gene entry with the largest delta, rounded to two decimals as spliceai prints."""
    gene, ds = max(per_gene, key=lambda g: max(g[1]))
    ds = tuple(round(v, 2) for v in ds)
    return gene, ds, max(ds)


def self_check(ann, models: list, records: list) -> None:
    """Batched scores must equal spliceai.utils.get_delta_scores (mask 0) to two decimals."""
    import spliceai.utils as su
    su.one_hot_encode = lambda s: one_hot(s).astype(np.int64)
    items = [(k, g, s, r, a) for k, rec in enumerate(records) for g, s, r, a in inputs(rec, ann)]
    batched = {}
    for key, gene, ds in score_batch(items, models):
        batched.setdefault(key, []).append((gene, ds))
    for k, rec in enumerate(records):
        ref = su.get_delta_scores(SimpleNamespace(chrom=rec.chrom, pos=rec.pos, ref=rec.ref, alts=[rec.alt]),
                                  ann, DIST, 0)
        expected = {e.split("|")[1]: tuple(float(v) for v in e.split("|")[2:6]) for e in ref}
        got = {g: tuple(round(v, 2) for v in ds) for g, ds in batched.get(k, [])}
        if any(abs(a - b) > 0.0101 for g in expected for a, b in zip(expected[g], got.get(g, (9, 9, 9, 9)))):
            raise SystemExit(f"self-check failed at {rec.chrom}:{rec.pos}: {expected} vs {got}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fasta", required=True, type=Path)
    ap.add_argument("--batch", type=int, default=32)
    a = ap.parse_args(argv)
    from spliceai.utils import Annotator
    ann = Annotator(str(a.fasta), "grch38")
    models = ann.models

    with gzip.open(VARIANTS, "rt", encoding="utf-8", newline="") as h:
        rows = [r for r in csv.DictReader(h, delimiter="\t")
                if r["new_label"] == "yes" and r["region"] == "deep_intronic" and r["variant_type"] == "snv"
                and r["stars"].isdigit() and int(r["stars"]) >= 1]
    study = {(r["chrom"], int(r["pos"])) for r in rows}
    check = []
    for k in range(40):                                   # positions inside BRCA1 introns, outside the study
        pos = 43_100_000 + 37 * k
        if ("17", pos) in study:
            continue
        base = ann.ref_fasta["chr17"][pos - 1].seq.upper()
        check.append(SimpleNamespace(chrom="17", pos=pos, ref=base, alt="A" if base != "A" else "G"))
    self_check(ann, models, check[:8])
    print("self-check passed on", len(check[:8]), "positions outside the study sets", flush=True)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    pending, results = [], {}
    recs = [SimpleNamespace(chrom=r["chrom"], pos=int(r["pos"]), ref=r["ref"], alt=r["alt"], vid=r["variation_id"])
            for r in rows]
    for n, rec in enumerate(recs, 1):
        pending += [(rec.vid, g, s, x, y) for g, s, x, y in inputs(rec, ann)]
        if len(pending) >= a.batch or n == len(recs):
            for vid, gene, ds in score_batch(pending, models) if pending else []:
                results.setdefault(vid, []).append((gene, ds))
            pending = []
            if n % 200 < 1 or n == len(recs):
                print(f"{n}/{len(recs)}", flush=True)
    with gzip.open(OUT, "wt", encoding="utf-8", newline="") as h:
        w = csv.writer(h, delimiter="\t", lineterminator="\n")
        w.writerow(("variation_id", "chrom", "pos", "ref", "alt", "spliceai_max", "DS_AG", "DS_AL", "DS_DG", "DS_DL",
                    "gene"))
        for rec in recs:
            if rec.vid in results:
                gene, ds, mx = best(results[rec.vid])
                w.writerow((rec.vid, rec.chrom, rec.pos, rec.ref, rec.alt, f"{mx:.2f}", *[f"{v:.2f}" for v in ds], gene))
            else:
                w.writerow((rec.vid, rec.chrom, rec.pos, rec.ref, rec.alt, "", "", "", "", "", ""))
    print("variants", len(recs), "scored", len(results))
    return 0


if __name__ == "__main__":
    sys.exit(main())
