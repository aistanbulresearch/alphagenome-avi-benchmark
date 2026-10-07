"""Part 2: experimental splicing measurements, downloaded and harmonised to GRCh38.

One function per dataset turns the authors' published table into rows of one common table
(`data/derived/part2/assay_variants.tsv.gz`). Every label is the authors' own call or the authors'
own threshold applied to the authors' own effect size; the continuous effect size is kept. Variants
are converted to GRCh38 (chrom without 'chr', 1-based pos, ref, alt on the forward strand; indels
left-aligned in VCF form) and every reference allele is checked against the GRCh38 FASTA; mismatches
are dropped and counted. hg19 positions are converted with `liftover.py`. No tool score is read.

Usage:
    python src/part2_data.py download      # fetch the public source files (download_public.py)
    python src/part2_data.py build         # write data/derived/part2/assay_variants.tsv.gz
    python src/part2_data.py build --fasta <hg38.fa> --summary <counts.json>
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import re
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from paths import DERIVED, EXPOSURE, RAW
from spreadsheet_readers import read_xlsx, parse_biff, read_xls, read_sheets, sheet_frame  # noqa: F401
from part2_genome import build_fai, Fasta, revcomp, normalize, variant_type, _chrom, _row, _label, _true  # noqa: F401

OUT = DERIVED / "part2" / "assay_variants.tsv.gz"
COLUMNS = ["dataset", "variant_id_in_source", "chrom", "pos", "ref", "alt", "variant_type", "effect_size",
           "effect_measure", "splice_altering", "threshold_rule", "gene", "region_note", "subset", "use"]

# Published path inside the raw folder -> public URL (authors' or journal's official copies).
_SPR18 = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-05748-7/MediaObjects/"
_SPR22 = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-31818-y/MediaObjects/"
_GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE307nnn/GSE307247/"
_COMPASS_GH = ("https://raw.githubusercontent.com/skoplik/ESL_MPRA_2025/"
               "67791aa433bf66ac0b7b095ede01d650c5de9cb1/Data_Pre-Processing/")
_MFASS_GH = ("https://raw.githubusercontent.com/KosuriLab/MFASS/"
             "9a8e4f27106be52aeb11acad27f95f5cded663a8/processed_data/snv/")
SOURCES = {
    "part2/kitzman/POU1F1_SpliceEffects.xls":
        "https://zenodo.org/api/records/7803918/files/POU1F1_SpliceEffects.xls/content",
    "part2/kitzman/WT1_SpliceEfects.xlsx":
        "https://zenodo.org/api/records/7803918/files/WT1_SpliceEfects.xlsx/content",
    "part2/kitzman/FormatManual.xlsx": "https://zenodo.org/api/records/7803918/files/FormatManual.xlsx/content",
    "part2/mst1r_braun2018/41467_2018_5748_MOESM3_ESM.pdf": _SPR18 + "41467_2018_5748_MOESM3_ESM.pdf",
    "part2/mst1r_braun2018/41467_2018_5748_MOESM6_ESM.xlsx": _SPR18 + "41467_2018_5748_MOESM6_ESM.xlsx",
    "part2/cd19_cortes_lopez2022/41467_2022_31818_MOESM3_ESM.pdf": _SPR22 + "41467_2022_31818_MOESM3_ESM.pdf",
    "part2/cd19_cortes_lopez2022/41467_2022_31818_MOESM5_ESM.xlsx": _SPR22 + "41467_2022_31818_MOESM5_ESM.xlsx",
    "part2/cd19_cortes_lopez2022/41467_2022_31818_MOESM6_ESM.xlsx": _SPR22 + "41467_2022_31818_MOESM6_ESM.xlsx",
    "part2/compass/GSE307247_column_descriptions_metadata.tsv.gz":
        _GEO + "suppl/GSE307247_column_descriptions_metadata.tsv.gz",
    "part2/compass/GSE307247_Processed_PSIs_All_Cells.csv.gz": _GEO + "suppl/GSE307247_Processed_PSIs_All_Cells.csv.gz",
    "part2/compass/GSE307247_family.soft.gz": _GEO + "soft/GSE307247_family.soft.gz",
    "part2/compass/github_67791aa/st_final_with_snp_and_coords_05_30_25.csv.gz":
        _COMPASS_GH + "st_final_with_snp_and_coords_05_30_25.csv.gz",
    "part2/opensplice/Supplementary_Table4_media-6.zip":
        "https://www.biorxiv.org/content/biorxiv/early/2026/05/23/2026.05.22.727141/DC6/embed/media-6.zip?download=true",
    "part2/mfass/github_9a8e4f2/snv_data_clean.txt": _MFASS_GH + "snv_data_clean.txt",
}
CLINVAR_VCF = "clinvar_20250615.vcf.gz"          # already in the raw folder (part 1 download)
_COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


# ---------------------------------------------------------------------------------------------------
# Datasets. Each harmonise_* takes the parsed source table and returns rows in source coordinates
# already converted to the forward strand; `finalize` checks them against the genome.
# ---------------------------------------------------------------------------------------------------
POU1F1_RULE = ("authors' calls: SDV = Z > 4.16 (Bonferroni p = 0.05) and >= 3-fold change from the null PSI "
               "of the beta, skip or other isoform in at least half the replicates; intermediate = Z met with a "
               "2- to 3-fold change, or SDV criteria met overall but in fewer than 7 replicates")
WT1_RULE = ("authors' calls: SDV = Bonferroni P < 0.05 and (SKIP or OTHER usage >= 20% above the null, or "
            "log2(KTS+/KTS-) >= 1.5 or <= -1); intermediate = same test and (>= 10% above the null, or log2 "
            "ratio >= 1 or <= -0.5)")
RON_RULE = ("authors' rule applied to the published model estimates: |change in isoform frequency| >= 5 points "
            "versus the wild-type mean and BH FDR < 0.05 (Stouffer-combined replicate z-scores, two-sided), "
            "for any of six isoforms")
RON_MEDIAN_RULE = ("no test possible (median-based estimate, not inferable by the model): |median change in "
                   "isoform frequency| >= 5 points for any isoform")
CD19_RULE = ("authors' call ('splicing-effective mutation' = Yes for any isoform): beyond the 2.5% and 97.5% "
             "quantiles of the wild-type minigene distribution in both replicates")
COMPASS_RULE = "authors' SDV rule: |delta logit(PSI)| >= 1 (mean of pooled-replicate values over cell lines)"
OPENSPLICE_RULE = "authors' call: Benjamini-Hochberg adjusted p < 0.1 (z-test on the logit scale, per exon)"
MFASS_RULE = ("authors' SDV call: delta inclusion index <= -0.5 for a wild-type exon with inclusion index "
              ">= 0.5 (column strong_lof)")
_HGVS = re.compile(r"^c\.([-*]?\d+)([+-]\d+)?([ACGT])>([ACGT])$")
POU1F1_COLUMNS = ["Chromosome", "Position (HG38)", "Position (HGVS)", "Ref", "Alt", "Variant Consequence",
                  "Splice Disruptive", "Intermediate", "Beta fold-change", "Skip fold-change", "Other fold-change"]
WT1_COLUMNS = ["HGVS Variant", "HG19 Position", "Variant Consequence", "Log2 KTS+/KTS- Ratio", "Splice Disruptive",
               "Intermediate"]


def harmonise_pou1f1(df: pd.DataFrame) -> list[dict]:
    """Gergics et al. 2021 Table S1. Alleles are given on the transcript (minus) strand."""
    out = []
    for r in df.to_dict("records"):
        if r.get("Position (HG38)") is None or r.get("Ref") is None:
            continue
        fc = [r.get(c) for c in ("Beta fold-change", "Skip fold-change", "Other fold-change")]
        fc = [float(x) for x in fc if isinstance(x, (int, float))]
        out.append(_row("kitzman_pou1f1_ex2", f"{r['Position (HGVS)']}{r['Ref']}>{r['Alt']}", r["Chromosome"],
                        int(r["Position (HG38)"]), revcomp(r["Ref"]), revcomp(r["Alt"]),
                        max(fc) if fc else np.nan,
                        "maximum fold change of beta, skip or other isoform PSI over the null (authors' metric)",
                        _label(_true(r["Splice Disruptive"]), _true(r["Intermediate"])), POU1F1_RULE, "POU1F1",
                        f"{r['Variant Consequence']}; {r['Position (HGVS)']}; NM_001122757.2",
                        "saturation_snv", "main"))
    return out


def harmonise_wt1(df: pd.DataFrame, lift) -> tuple[list[dict], int]:
    """Smith et al. 2023 Supplementary Table 2: HGVS on the transcript (minus strand) and hg19 position."""
    out, unlifted = [], 0
    for r in df.to_dict("records"):
        m = _HGVS.match(str(r.get("HGVS Variant") or "").strip())
        if not m or r.get("HG19 Position") is None:
            continue
        hits = lift(("11", int(r["HG19 Position"]), int(r["HG19 Position"])))
        if len(hits) != 1:
            unlifted += 1
            continue
        ratio = r.get("Log2 KTS+/KTS- Ratio")
        out.append(_row("kitzman_wt1_ex9", r["HGVS Variant"], hits[0][0], hits[0][1], revcomp(m.group(3)),
                        revcomp(m.group(4)), float(ratio) if isinstance(ratio, (int, float)) else np.nan,
                        "log2(KTS+/KTS-) isoform ratio, signed (SKIP and OTHER usage not included)",
                        _label(_true(r["Splice Disruptive"]), _true(r["Intermediate"])), WT1_RULE, "WT1",
                        f"{r['Variant Consequence']}; {r['HGVS Variant']}", "saturation_snv", "main"))
    return out, unlifted


_MUT = re.compile(r"^([ACGT]+)(\d+)([ACGT]+)\*?$")
RON_ISOFORMS = ["AE inclusion (%)", "AE skipping (%)", "full IR (%)", "first IR (%)", "second IR (%)", "other (%)"]


def bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out


def ron_calls(reps: list[pd.DataFrame]) -> pd.DataFrame:
    """Per mutation: max |change| (points) and the authors' call, from replicate sheets.

    Each replicate frame has the rows 'wt mean', 'wt SD' (any case) and one row per mutation, indexed
    by mutation, with the six isoform columns. z = (estimate - WT mean) / WT SD per replicate,
    Stouffer-combined, two-sided p, BH within each isoform; change = mean over replicates of
    (estimate - WT mean). Only mutations present in every replicate are tested.
    """
    norm = []
    for df in reps:
        low = {str(i).strip().lower(): i for i in df.index}
        mean = df.loc[low["wt mean"], RON_ISOFORMS].astype(float)
        sd = df.loc[low["wt sd"], RON_ISOFORMS].astype(float)
        body = df.drop(index=[low["wt mean"], low["wt sd"]])[RON_ISOFORMS].astype(float)
        norm.append((body, mean, sd))
    common = norm[0][0].index
    for body, _, _ in norm[1:]:
        common = common.intersection(body.index)
    z = sum((b.loc[common] - m) / s for b, m, s in norm) / math.sqrt(len(norm))
    delta = sum(b.loc[common] - m for b, m, _ in norm) / len(norm)
    p = z.apply(lambda col: pd.Series([math.erfc(abs(v) / math.sqrt(2)) for v in col], index=col.index))
    q = p.apply(lambda col: pd.Series(bh(col.values), index=col.index))
    sig = ((q < 0.05) & (delta.abs() >= 5)).any(axis=1)
    return pd.DataFrame({"max_abs_change": delta.abs().max(axis=1), "splice_altering": sig})


def _minigene_to_genome(note: str, offset: int, strand: str, first: int | None = None
                        ) -> tuple[int, str, str] | None:
    """'A26T' / 'CT38C' (minigene coordinates, transcript orientation) -> (start, ref, alt) on the
    forward strand of the source build; strand '+' maps p to offset + p, '-' maps p to offset - p.
    `first` is the minigene position of the first base of the ref string when it differs from the
    number in the notation."""
    m = _MUT.match(note)
    if not m:
        return None
    ref, alt = m.group(1), m.group(3)
    p = int(m.group(2)) if first is None else first
    if strand == "+":
        return offset + p, ref, alt
    return offset - (p + len(ref) - 1), revcomp(ref), revcomp(alt)


def harmonise_ron(reps: list[pd.DataFrame], average: pd.DataFrame, median: pd.DataFrame, lift,
                  fetch19=None) -> tuple[list[dict], dict]:
    """Braun et al. 2018 Supplementary Data 3, HEK293T. hg19, minus strand; hg19 = C - minigene position.

    The constant C is taken from the SNV rows (genomic position + minigene position) and must be
    unique. The authors' genomic position is the first base of the ref string (the anchor base of an
    indel); the number in the indel notation is the anchor in one sheet and the deleted base in the
    other, so the anchor is taken from the genomic position (C - position). Positions are lifted to
    GRCh38 by the offset of the lifted first base.
    """
    snv = average[average["mutation type"] == "SNV"]
    consts = {int(str(g).split(":")[1]) + int(_MUT.match(m).group(2))
              for m, g in zip(snv.index, snv["genomic position"])}
    if len(consts) != 1:
        raise ValueError(f"RON minigene-to-hg19 offset not unique: {sorted(consts)[:5]}")
    const = consts.pop()
    calls = ron_calls(reps)
    out, info = [], {"hg19_offset": const, "not_tested_in_all_replicates": 0, "unlifted": 0}

    def add(note, mtype, genomic, effect, label, rule, subset):
        g = _minigene_to_genome(note, const, "-", first=const - int(str(genomic).split(":")[1]))
        if g is None:
            return
        start, ref, alt = g
        hits = lift(("3", start, start))
        if len(hits) != 1:
            info["unlifted"] += 1
            return
        out.append(_row("braun2018_mst1r_ex11", note, hits[0][0], hits[0][1], ref, alt, effect,
                        "maximum |change in isoform frequency| over six isoforms, percentage points "
                        "(model-inferred single-mutation effect)", _label(label), rule, "MST1R",
                        f"minigene position {_MUT.match(note).group(2)} (exons 10-12 minigene); type {mtype}",
                        subset, "main"))

    for note, r in average.iterrows():
        if r["mutation type"] not in ("SNV", "DEL", "IN", "INS"):
            continue
        if note not in calls.index:
            info["not_tested_in_all_replicates"] += 1
            continue
        c = calls.loc[note]
        add(note, r["mutation type"], r["genomic position"], float(c["max_abs_change"]),
            bool(c["splice_altering"]), RON_RULE, "model_inferred")
    cols = [c for c in median.columns if str(c).startswith("Median delta")]
    for _, r in median.iterrows():
        eff = float(np.nanmax(np.abs(r[cols].astype(float).values)))
        add(str(r["Mutation"]), r["Mutation type"], r["Genomic position"], eff, eff >= 5, RON_MEDIAN_RULE,
            "median_estimate")
    return out, info


def cd19_offset(validated: pd.DataFrame) -> int:
    """hg38 position minus minigene position, from the authors' validated mutations (must be unique)."""
    consts = set()
    for _, r in validated.iterrows():
        coord, note = r.get("coordinate (hg38)"), r.get("corresponding minigene mutation")
        if isinstance(coord, str) and isinstance(note, str) and _MUT.match(note.strip()):
            consts.add(int(coord.split(":")[1]) - int(_MUT.match(note.strip()).group(2)))
    if len(consts) != 1:
        raise ValueError(f"CD19 minigene offset not unique: {consts}")
    return consts.pop()


def harmonise_cd19(effects: pd.DataFrame, offset: int) -> list[dict]:
    """Cortes-Lopez et al. 2022 Supplementary Data 3 ('Mutation effects'): six rows per mutation."""
    out = []
    for note, g in effects.groupby("mutation", sort=False):
        pos = _minigene_to_genome(str(note), offset, "+")
        if pos is None:
            continue
        eff = pd.to_numeric(g["average delta frequency"], errors="coerce").abs().max()
        out.append(_row("cortes_lopez2022_cd19", note, "16", pos[0], pos[1], pos[2], float(eff),
                        "maximum |average change in isoform frequency| over six isoform classes, fraction "
                        "(model-inferred single-mutation effect)",
                        _label((g["splicing-effective mutation"].astype(str).str.strip() == "Yes").any()),
                        CD19_RULE, "CD19", f"{g['region'].iloc[0]}; minigene position {_MUT.match(str(note)).group(2)}",
                        "model_inferred", "main"))
    return out


COMPASS_CELLS = ("HEK", "HeLa", "K562", "MCF7", "HMC3")


def harmonise_compass(df: pd.DataFrame, clinvar_keys: set) -> tuple[list[dict], dict]:
    """Koplik et al. 2025, GEO processed table. Single variants only; doubles and references are counted.

    subset: exac / geuvadis (authors' source labels); the authors' 'clinvar_single' source mixes ClinVar
    variants and random SNVs, so a 'clinvar_single' allele present in the given ClinVar release is the
    ClinVar subset and the rest are the random SNVs. Effect = mean over cell lines of the pooled
    delta logit; a variant whose reference PSI is exactly 0 or 1 in every cell line where it was
    measured is not used (the authors exclude such sequences when evaluating models).
    """
    info = {"references": int((df["snp"] == "none").sum()), "doubles": int((df["seq_type"] == "double").sum()),
            "single_not_measured": 0}
    singles = df[(df["seq_type"] == "single") & df["variant_hg38"].notna()]
    dl = singles[[f"{c}_delta_logit_pooled" for c in COMPASS_CELLS]].apply(pd.to_numeric, errors="coerce")
    ref = singles[[f"{c}_wt_pooled_psi_raw" for c in COMPASS_CELLS]].apply(pd.to_numeric, errors="coerce")
    measured = dl.notna().values
    extreme = ((ref.values == 0) | (ref.values == 1)) | ~measured
    mean = dl.mean(axis=1, skipna=True)
    out = []
    for i, r in enumerate(singles.to_dict("records")):
        if not measured[i].any():
            info["single_not_measured"] += 1
            continue
        chrom, pos, change = str(r["variant_hg38"]).split(":")
        ref_a, alt_a = change.split(">")
        key = (_chrom(chrom), int(pos), ref_a, alt_a)
        if r["source"] == "exac":
            subset, use = "exac", "main"
        elif r["source"] == "geuvadis":
            subset, use = "geuvadis", "not_used"
        elif key in clinvar_keys:
            subset, use = "clinvar", "not_used"
        else:
            subset, use = "random_snv", "main"
        if use == "main" and extreme[i].all():
            subset, use = subset + ";reference_psi_0_or_1", "not_used"
        start, end = int(r["exon_start_hg38"]), int(r["exon_end_hg38"])
        where = "exon" if start <= int(pos) <= end else "intron"
        cells = ",".join(c for c, ok in zip(COMPASS_CELLS, measured[i]) if ok)
        gene = str(r["gene_exon"]).split(" exon")[0]
        out.append(_row("compass", f"{r['Reference']}:{r['snp']}", chrom, pos, ref_a, alt_a, float(mean.iloc[i]),
                        "delta logit(PSI), mean over cell lines of pooled-replicate values (cells: " + cells + ")",
                        _label(abs(mean.iloc[i]) >= 1), COMPASS_RULE, gene,
                        f"{where}; {r['gene_exon']}; {r['transcript_id']}", subset, use))
    return out, info


def harmonise_opensplice(df: pd.DataFrame) -> tuple[list[dict], dict]:
    """Quarantani et al. 2026 Supplementary Table 4 (hg38; deletions in VCF form with a 'del' position)."""
    info = {"wt_rows": int((df["mut_type"] == "wt").sum()), "not_measured": 0}
    out = []
    for r in df[df["mut_type"] != "wt"].to_dict("records"):
        d = r.get("delta_logit")
        if d is None or (isinstance(d, float) and math.isnan(d)) or pd.isna(r.get("psi")):
            info["not_measured"] += 1
            continue
        pos_s = str(r["POS"])
        if pos_s.endswith("del"):
            pos = int(pos_s[:-3].split("_")[0]) - 1          # anchor base precedes the deleted span
        else:
            pos = int(float(pos_s))
        gene = str(r["exon_id"]).rsplit("_e", 1)[0]
        out.append(_row("opensplice", r["variant_id"], r["CHROM"], pos, r["REF"], r["ALT"], float(d),
                        "delta logit(PSI) (mutant - WT), minigene",
                        _label(str(r["significant"]).strip().lower() == "yes"), OPENSPLICE_RULE, gene,
                        f"{r['region']}; {r['exon_id']}; library {r['lib_id']}", "saturation", "main"))
    return out, info


def harmonise_mfass(df: pd.DataFrame) -> tuple[list[dict], dict]:
    """Cheung et al. 2019, snv_data_clean.txt: ExAC SNVs ('mutant'), hg38 forward-strand alleles."""
    mut = df[df["category"] == "mutant"]
    called = mut[mut["strong_lof"].notna()]
    info = {"mutant_rows": len(mut), "no_call": len(mut) - len(called)}
    out = []
    for r in called.to_dict("records"):
        out.append(_row("mfass", r["id"], r["chr"], int(r["snp_position_hg38_1based"]), r["ref_allele"],
                        r["alt_allele"], float(r["v2_dpsi"]), "delta exon inclusion index (library v2)",
                        _label(_true(r["strong_lof"])), MFASS_RULE, r["ensembl_id"],
                        f"{r['label']}; exon {r['ensembl_id']}; tested by the developers of compared tools",
                        "exac", "comparison"))
    return out, info


# ---------------------------------------------------------------------------------------------------
# Checks, de-duplication and assembly
# ---------------------------------------------------------------------------------------------------
def finalize(rows: list[dict], fetch) -> tuple[pd.DataFrame, dict]:
    """Check reference alleles, left-align indels, flag duplicate alleles within a dataset.

    A row whose reference allele differs from the genome is dropped and counted. When the same
    allele occurs more than once within a dataset among usable rows, the first row in
    (variant_id_in_source) order keeps its `use`; the others become not_used with ';duplicate'.
    """
    kept, mismatches = [], {}
    for r in rows:
        ref, alt = r["ref"], r["alt"]
        if fetch(r["chrom"], r["pos"], r["pos"] + len(ref) - 1) != ref:
            mismatches[r["dataset"]] = mismatches.get(r["dataset"], 0) + 1
            continue
        if not (len(ref) == 1 and len(alt) == 1):
            r["pos"], r["ref"], r["alt"] = normalize(fetch, r["chrom"], r["pos"], ref, alt)
        r["variant_type"] = variant_type(r["ref"], r["alt"])
        kept.append(r)
    df = pd.DataFrame(kept, columns=COLUMNS)
    if df.empty:
        return df, {"ref_mismatch": mismatches, "duplicates": {}}
    df = df.sort_values(["dataset", "variant_id_in_source"], kind="stable").reset_index(drop=True)
    usable = df["use"] != "not_used"
    dup = usable & df[usable].duplicated(["dataset", "chrom", "pos", "ref", "alt"]).reindex(df.index,
                                                                                            fill_value=False)
    df.loc[dup, "subset"] = df.loc[dup, "subset"] + ";duplicate"
    df.loc[dup, "use"] = "not_used"
    dups = df[dup].groupby("dataset").size().to_dict()
    df = df.sort_values(["dataset", "chrom", "pos", "ref", "alt", "variant_id_in_source"], kind="stable")
    return df.reset_index(drop=True), {"ref_mismatch": mismatches, "duplicates": dups}


def clinvar_keys(path: Path, wanted: set) -> set:
    """(chrom, pos, ref, alt) of ClinVar VCF records at the wanted (chrom, pos) positions."""
    keys = set()
    with gzip.open(path, "rt", encoding="utf-8") as h:
        for line in h:
            if line.startswith("#"):
                continue
            c, p, _, ref, alt = line.split("\t", 5)[:5]
            if (c, int(p)) in wanted:
                for a in alt.split(","):
                    keys.add((c, int(p), ref, a))
    return keys


def load_all(raw: Path, lift, clinvar_vcf: Path | None) -> tuple[list[dict], dict]:
    """Read every source file and return harmonised rows (before the genome check) and counts."""
    rows, info = [], {}
    kz = raw / "part2" / "kitzman"
    pou = sheet_frame(read_sheets(kz / "POU1F1_SpliceEffects.xls")["supp_tbl.2021-0420"], 2, POU1F1_COLUMNS)
    r = harmonise_pou1f1(pou)
    rows += r
    info["kitzman_pou1f1_ex2"] = {"source_rows": len(pou), "harmonised": len(r)}
    wt1 = sheet_frame(read_sheets(kz / "WT1_SpliceEfects.xlsx")["Supplementary Table 2"], 1, WT1_COLUMNS)
    r, unlifted = harmonise_wt1(wt1, lift)
    rows += r
    info["kitzman_wt1_ex9"] = {"source_rows": len(wt1), "harmonised": len(r), "unlifted": unlifted}

    ron = read_sheets(raw / "part2" / "mst1r_braun2018" / "41467_2018_5748_MOESM6_ESM.xlsx")

    def ron_sheet(name: str, index: str) -> pd.DataFrame:
        return sheet_frame(ron[name], 0).set_index(index)

    reps = [ron_sheet(f"HEK293T rep{i}", "mutation") for i in (1, 2, 3)]
    avg = ron_sheet("HEK293T average frequency", "mutation")
    avg = avg.drop(index=[i for i in avg.index if str(i).lower().startswith("wt ")])
    med = sheet_frame(ron["HEK293T median-based estimates"], 0)
    r, ron_info = harmonise_ron(reps, avg, med, lift)
    rows += r
    info["braun2018_mst1r_ex11"] = {"source_rows": len(avg) + len(med), "harmonised": len(r), **ron_info}

    cd = read_sheets(raw / "part2" / "cd19_cortes_lopez2022" / "41467_2022_31818_MOESM6_ESM.xlsx")
    effects = sheet_frame(cd["Mutation effects"], 0)
    offset = cd19_offset(sheet_frame(cd["19 validated mutations"], 0))
    r = harmonise_cd19(effects, offset)
    rows += r
    info["cortes_lopez2022_cd19"] = {"source_rows": int(effects["mutation"].nunique()), "harmonised": len(r),
                                     "hg38_offset": offset}

    comp = pd.read_csv(raw / "part2" / "compass" / "GSE307247_Processed_PSIs_All_Cells.csv.gz", low_memory=False,
                       usecols=lambda c: c in ("Reference", "snp", "source", "seq_type", "gene_exon", "transcript_id",
                                               "exon_start_hg38", "exon_end_hg38", "variant_hg38")
                       or c.endswith(("_delta_logit_pooled", "_wt_pooled_psi_raw")))
    singles = comp[(comp["seq_type"] == "single") & (comp["source"] == "clinvar_single")]
    wanted = {(f"{c[3:] if c.startswith('chr') else c}", int(p)) for c, p, _ in
              (str(v).split(":") for v in singles["variant_hg38"].dropna())}
    keys = clinvar_keys(clinvar_vcf, wanted) if clinvar_vcf else set()
    r, comp_info = harmonise_compass(comp, keys)
    rows += r
    info["compass"] = {"source_rows": len(comp), "harmonised": len(r), **comp_info,
                       "clinvar_release_for_subset": clinvar_vcf.name if clinvar_vcf else None}

    with zipfile.ZipFile(raw / "part2" / "opensplice" / "Supplementary_Table4_media-6.zip") as z:
        osp = pd.read_csv(z.open("Supplementary_Table4.tsv"), sep="\t", low_memory=False,
                          usecols=["CHROM", "POS", "REF", "ALT", "exon_id", "variant_id", "lib_id", "mut_type",
                                   "region", "psi", "delta_logit", "significant"])
    r, osp_info = harmonise_opensplice(osp)
    rows += r
    info["opensplice"] = {"source_rows": len(osp), "harmonised": len(r), **osp_info}

    mf = pd.read_csv(raw / "part2" / "mfass" / "github_9a8e4f2" / "snv_data_clean.txt", sep="\t", low_memory=False,
                     usecols=["id", "ensembl_id", "chr", "ref_allele", "alt_allele", "snp_position_hg38_1based",
                              "label", "category", "v2_dpsi", "strong_lof"])
    r, mf_info = harmonise_mfass(mf)
    rows += r
    info["mfass"] = {"source_rows": len(mf), "harmonised": len(r), **mf_info}
    return rows, info


def summarise(df: pd.DataFrame) -> dict:
    out = {}
    for ds, g in df.groupby("dataset"):
        d = {"rows": len(g), "by_use": g["use"].value_counts().to_dict(),
             "by_type": g["variant_type"].value_counts().to_dict()}
        for use, gu in g.groupby("use"):
            d[f"labels_{use}_SNV"] = gu[gu["variant_type"] == "SNV"]["splice_altering"].value_counts().to_dict()
            d[f"labels_{use}_indel"] = gu[gu["variant_type"] != "SNV"]["splice_altering"].value_counts().to_dict()
        d["by_subset"] = g["subset"].value_counts().to_dict()
        out[ds] = d
    return out


def write(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    out = df.copy()
    out["effect_size"] = out["effect_size"].map(lambda v: "" if pd.isna(v) else f"{v:.6g}")
    with gzip.GzipFile(path, "wb", mtime=0) as gz:
        gz.write(out[COLUMNS].to_csv(sep="\t", index=False, lineterminator="\n").encode("utf-8"))


def download_all(raw: Path) -> None:
    import download_public
    for name, url in SOURCES.items():
        if (raw / name).exists():
            print(f"present: {name}")
            continue
        download_public.main(["--name", name, url])


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("command", choices=["download", "build"])
    ap.add_argument("--fasta", default=str(RAW / "hg38.fa"))
    ap.add_argument("--chain", default=str(EXPOSURE / "hg19ToHg38.over.chain.gz"))
    ap.add_argument("--clinvar", default=str(RAW / CLINVAR_VCF))
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--summary", help="optional JSON file for the counts")
    args = ap.parse_args(argv)
    if args.command == "download":
        download_all(RAW)
        return 0
    import liftover
    fasta = Fasta(args.fasta)
    rows, info = load_all(RAW, liftover.load(Path(args.chain)), Path(args.clinvar))
    df, checks = finalize(rows, fasta.fetch)
    write(df, Path(args.out))
    report = {"inputs": info, "checks": checks, "table": summarise(df)}
    text = json.dumps(report, indent=2, default=str)
    if args.summary:
        Path(args.summary).write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0


__all__ = ["COLUMNS", "SOURCES", "Fasta", "build_fai", "normalize", "variant_type", "revcomp", "read_xlsx",
           "read_xls", "parse_biff", "read_sheets", "sheet_frame", "bh", "ron_calls", "harmonise_pou1f1",
           "harmonise_wt1", "harmonise_ron", "cd19_offset", "harmonise_cd19", "harmonise_compass",
           "harmonise_opensplice", "harmonise_mfass", "finalize", "clinvar_keys", "load_all", "main"]

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
