"""One table for the clinical-label analysis: variant list, gnomAD frequencies, submission flags and the
tool scores, joined by VariationID; plus the set and region definitions of STUDY_DESIGN.md section 4.4.
"""
from __future__ import annotations

import csv
import gzip
from pathlib import Path

import numpy as np
import pandas as pd

from paths import DERIVED

VARIANTS = DERIVED / "variants"
SCORES = DERIVED / "scores"
RARE = 0.001
MAIN_REGION = ("donor_region", "acceptor_region", "near_intronic", "deep_intronic")
TOOLS = {                      # name: (file, column); higher score = more pathogenic
    "AVI": ("avi.tsv.gz", "avi_raw"),
    "SpliceAI": ("spliceai.tsv.gz", "spliceai_max"),
    "SpliceAI masked": ("spliceai_masked.tsv.gz", "spliceai_max"),
    "SpliceAI 4999": ("spliceai_d4999.tsv.gz", "spliceai_max"),
    "CADD": ("cadd.tsv.gz", "cadd_phred"),
    "phyloP 100-way": ("conservation.tsv.gz", "phyloP100"),
    "phastCons 470-way": ("conservation.tsv.gz", "phastCons470"),
    "phyloP 241-way": ("conservation.tsv.gz", "phyloP241"),
}


def _read(path: Path, usecols=None) -> pd.DataFrame:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as h:
        return pd.read_csv(h, sep="\t", dtype=str, usecols=usecols, keep_default_na=False)


def rare_flags(g: pd.DataFrame) -> pd.DataFrame:
    """rare: joint FAF (group maximum, 95%) below 0.001, else joint AF below 0.001; absent = rare."""
    faf = pd.to_numeric(g["fafmax_faf95_max_joint"], errors="coerce")
    af = pd.to_numeric(g["AF_joint"], errors="coerce")
    present = g["in_gnomad"] == "yes"
    freq = faf.where(faf.notna(), af)
    rare = ~present | freq.isna() | (freq < RARE)
    return pd.DataFrame({"variation_id": g["variation_id"], "in_gnomad": present, "rare": rare})


def load() -> pd.DataFrame:
    calib = set(pd.read_csv(VARIANTS / "calibration_sample.tsv", sep="\t", dtype=str)["variation_id"])
    cols = ["variation_id", "chrom", "pos", "ref", "alt", "gene_id", "label", "stars", "chromosome_group",
            "variant_type", "family", "region", "intron_offset", "mane_select", "new_label",
            "inside_long_unknown_record"]
    v = _read(VARIANTS / "clinvar_variants.tsv.gz", cols)
    v = v[(v["variant_type"] == "snv") & ((v["new_label"] == "yes") | v["variation_id"].isin(calib))]
    v = v[v["family"].isin(["intronic", "utr"]) | (v["region"] == "synonymous")]
    gn = pd.concat([_read(VARIANTS / "gnomad_v4.1_joint.tsv.gz"),
                    _read(VARIANTS / "gnomad_v4.1_joint_calibration.tsv.gz")]).drop_duplicates("variation_id")
    v = v.merge(rare_flags(gn), on="variation_id", how="left")
    fl = _read(VARIANTS / "submission_flags.tsv.gz")
    v = v.merge(fl[["variation_id", "mentions_spliceai", "mentions_any_splice_predictor"]], on="variation_id",
                how="left")
    for tool, (name, col) in TOOLS.items():
        path = SCORES / name
        if path.exists():
            s = _read(path, ["variation_id", col]).rename(columns={col: tool})
            s[tool] = pd.to_numeric(s[tool], errors="coerce")
            v = v.merge(s, on="variation_id", how="left")
    off = pd.to_numeric(v["intron_offset"], errors="coerce")
    v["Distance"] = -off.abs()                     # nearer the exon boundary = more likely pathogenic
    v["y"] = (v["label"] == "P/LP").astype(int)
    v["stars_n"] = pd.to_numeric(v["stars"], errors="coerce")
    v["calibration"] = v["variation_id"].isin(calib)
    return v.reset_index(drop=True)


def tools_present(v: pd.DataFrame) -> list:
    return [t for t in list(TOOLS) + ["Distance"] if t in v.columns]


def sets(v: pd.DataFrame) -> dict:
    """Boolean masks of the analysis sets (section 4.4)."""
    new = (v["new_label"] == "yes") & (v["stars_n"] >= 1) & (v["mane_select"] == "yes")
    uv = v["chromosome_group"].isin(["untouched", "validation"])
    rare = v["rare"].fillna(True).astype(bool)
    main = new & uv & rare
    return {
        "main": main,
        "all frequencies": new & uv,
        "present and rare": main & v["in_gnomad"].fillna(False).astype(bool),
        "training chromosomes": new & (v["chromosome_group"] == "training") & rare,
        "two stars": main & (v["stars_n"] >= 2),
        "stricter exclusion": main & (v["inside_long_unknown_record"] == "no"),
        "predictor mentioned": main & (v["mentions_any_splice_predictor"] == "yes"),
        "predictor not mentioned": main & (v["mentions_any_splice_predictor"] == "no"),
        "calibration sample": v["calibration"] & rare & (v["stars_n"] >= 1),
    }


def regions(v: pd.DataFrame) -> dict:
    r = v["region"]
    return {
        "intronic >=3 bp (main)": r.isin(MAIN_REGION),
        "donor region": r == "donor_region",
        "acceptor region": r == "acceptor_region",
        "near intronic": r == "near_intronic",
        "deep intronic": r == "deep_intronic",
        "canonical": r == "canonical",
        "UTR": r.isin(["utr5", "utr3"]),
        "5' UTR": r == "utr5",
        "3' UTR": r == "utr3",
        "synonymous": r == "synonymous",
    }


__all__ = ["load", "sets", "regions", "tools_present", "MAIN_REGION", "TOOLS", "rare_flags"]
