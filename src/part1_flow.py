"""Variant flow of Part 1, from the ClinVar release to the main set, with P/LP and B/LB counts per step.
Filters from the new-label set to the main set apply in the order listed; each row counts the variants that
step removes. Writes results/part1/flow.tsv (used by the study design figure and a supplementary table).

Usage:
    python src/part1_flow.py
"""
from __future__ import annotations

import sys

import pandas as pd

import part1_data as pdata
from paths import DERIVED, ROOT

OUT = ROOT / "results" / "part1" / "flow.tsv"
EXCLUSIONS = [("labelled before", "already P/LP or B/LB at or before the 2025-06-15 release (or the 2024-07, 2025-07 archives)"),
              ("developers' benchmark or submitter", "at AVI developers' benchmark positions or from a development submitter"),
              ("unresolved historical identity (same position)", "older record of unresolved identity at the same position"),
              ("identity collision", "identity collisions")]


def counts(labels: pd.Series) -> dict:
    return {"P/LP": int((labels == "P/LP").sum()), "B/LB": int((labels == "B/LB").sum()), "total": int(len(labels))}


def flow() -> pd.DataFrame:
    d = pd.read_csv(DERIVED / "variants" / "clinvar_variants.tsv.gz", sep="\t", dtype=str, keep_default_na=False,
                    usecols=["variation_id", "label", "new_label", "exclusion", "variant_type", "family", "region"])
    rows = [{"stage": "release", "step": "ClinVar October 2026, P/LP or B/LB", **counts(d["label"])}]
    for key, text in EXCLUSIONS:
        hit = d["exclusion"].str.startswith(key)
        rows.append({"stage": "excluded", "step": text, **counts(d.loc[hit, "label"])})
    new = d[d["new_label"] == "yes"]
    rows.append({"stage": "new-label set", "step": "new-label set", **counts(new["label"])})
    snv = new[new["variant_type"] == "snv"]
    rows.append({"stage": "not in the main set", "step": "not SNV", **counts(new.loc[new["variant_type"] != "snv", "label"])})
    in_regions = snv["family"].isin(["intronic", "utr"]) | (snv["region"] == "synonymous")
    unparsed = snv["region"] == "unparsed"
    rows.append({"stage": "not in the main set", "step": "coding other than synonymous (missense, nonsense, other)",
                 **counts(snv.loc[~in_regions & ~unparsed, "label"])})
    rows.append({"stage": "not in the main set", "step": "HGVS name not parsed",
                 **counts(snv.loc[~in_regions & unparsed, "label"])})
    v = pdata.load()
    v = v[v["new_label"] == "yes"]
    if len(v) != int(in_regions.sum()):
        raise ValueError("the loaded rows do not match the study-region SNVs of the new-label set")
    keep = pd.Series(True, index=v.index)
    for step, ok, note in (("no review star", v["stars_n"] >= 1, ""),
                           ("no MANE Select name", v["mane_select"] == "yes", ""),
                           ("on training chromosomes", v["chromosome_group"].isin(["untouched", "validation"]),
                            "supporting analysis"),
                           ("not rare in gnomAD v4.1 (FAF95 >= 0.001)", v["rare"].fillna(True).astype(bool),
                            "supporting analysis")):
        removed = keep & ~ok
        rows.append({"stage": "not in the main set", "step": step, "note": note, **counts(v.loc[removed, "label"])})
        keep &= ok
    rows.append({"stage": "main set", "step": "main set", **counts(v.loc[keep, "label"])})
    return pd.DataFrame(rows).fillna({"note": ""})


def main(argv: list[str] | None = None) -> int:
    out = flow()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False)
    print(out.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
