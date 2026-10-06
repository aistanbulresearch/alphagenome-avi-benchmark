# Independent evaluation of the AlphaGenome Variant Impact score

Study design, 2026-10-07.

## 1. Question

How well does the AlphaGenome Variant Impact (AVI) score identify disease-causing non-coding and splice-affecting variants on data its developers did not use, compared with established tools? And is it calibrated well enough for use as clinical evidence?

We answer this with two kinds of ground truth:

1. **Clinical labels.** ClinVar pathogenic and benign classifications added after the developers' ClinVar snapshot.
2. **Experimental measurements.** Splicing assays that tested variants without selecting them by any tool's prediction.

A third part, the prioritization of causal variants in real genomes, is optional. It is decided after parts 1 and 2 (section 7).

## 2. Integrity principle

There is one rule: **no choice about data or analysis is made after looking at tool scores.**

In practice:

- The variant sets, region definitions and metrics in this document are written down and committed before any score is retrieved. The commit date is the record.
- Every metric we compute is reported, whether it favours AVI or not.
- Code, variant lists and score files (where licences allow) are released, so that anyone can reproduce every number.

Counts of variants (how many there are per region) contain no scores. They may be used to decide which analyses are reported as main and which as supporting.

## 3. Evaluated tools

| Role | Tool |
|---|---|
| Evaluated | AVI (AlphaGenome Variant Impact score) |
| Splicing comparators | SpliceAI, Pangolin |
| Genome-wide comparators | CADD v1.7 |
| Conservation | phyloP, PhastCons |
| Simple baselines | distance to the nearest splice site |

Before scoring, we read each tool's documentation and terms of use. We note which datasets each tool was trained or evaluated on. A dataset that a compared tool was trained on is not used to compare that tool.

AVI uses conservation scores and AlphaMissense as inputs. The conservation comparison is therefore reported as "what AVI adds beyond conservation", not as a head-to-head comparison.

## 4. Part 1: clinical labels

### 4.1 Variant set

- **Source:** ClinVar monthly releases. The developers used the ClinVar release of 2025-06-15.
- **Main set ("new labels"):**
  - The variant has a pathogenic/likely pathogenic or benign/likely benign classification in the October 2026 release.
  - It had no such classification in the 2025-06-15 release or in the earlier archives used for history (2024-07, 2025-07).
  - It has at least one review star.
- **Excluded from the main set:**
  - Variants at positions in benchmark datasets the developers used for training, selection or evaluation (saturation genome editing screens, TraitGym, ncVarDB).
  - Submissions from the organisation that supplied the developers' GREGoR cases.
- **Comparison set ("all labels"):** all ClinVar-labelled variants of the same regions. The developers evaluated AVI on these labels. Results on this set are reported next to the main set and labelled accordingly.
- **Chromosome groups:** results are broken down by the developers' training, validation and untouched chromosome groups.
- **Rare variants:** variants that are common in gnomAD v4.1 are reported separately, because AVI was trained on population frequency.

The variant counts of this set were computed from the October 2026 release before any score was retrieved.

### 4.2 Regions

The regions are assigned from the HGVS description on the MANE Select transcript:

- canonical splice site (±1, ±2);
- splice region;
- near intronic;
- deep intronic;
- 5' UTR;
- 3' UTR;
- synonymous;
- other coding.

Main reporting is by region. Neighbouring regions are pooled when a region is too small to estimate a metric, for example intronic ≥3 bp from the exon, or UTR pooled.

### 4.3 Metrics

- **Discrimination:** AUROC and AUPRC per region for every tool, and the paired AUROC difference between AVI and SpliceAI (intronic regions) or CADD (UTR). 95% intervals come from a bootstrap that resamples genes, because variants of the same gene are not independent.
- **Clinical evidence strength:** we use ClinGen's calibration framework (Pejaver et al. 2022). For each tool we find the score thresholds at which the positive likelihood ratio reaches supporting (2.406), moderate (5.790) and strong (33.53) evidence, with confidence intervals. Thresholds are derived on the comparison set and checked on the main set, so the main set is not used to choose them.
- **Coverage:** the fraction of variants each tool can score.

Small regions are reported with their intervals. Wide intervals are a result, not a reason to drop a region.

## 5. Part 2: experimental splicing measurements

### 5.1 Data

We use selection-free assays: libraries whose variants were chosen without reference to any tool's prediction. Candidates still need checking against their publications:

- saturation minigene assays: Kitzman lab POU1F1 and WT1; MST1R/RON (Braun et al. 2018); CD19 (Cortés-López et al. 2022);
- COMPASS (GEO GSE307247, random and ExAC subsets);
- OpenSplice (bioRxiv 2026);
- MFASS. AlphaGenome and Pangolin report results on it, so it is a comparison dataset, not a main one.

Each dataset's own threshold for "splice-altering" is taken from its publication. A variant present in both parts is counted in one part only.

### 5.2 Metrics

- AUROC for splice-altering versus not, per dataset and pooled, with dataset-level intervals.
- Rank correlation with the measured effect size.

## 6. Repository contents

The repository holds research material only:

```
README.md          what the study is and how to reproduce it
STUDY_DESIGN.md    this document
src/               code
tests/             tests of the code
data/              download list with sources, dates and checksums; derived variant lists (no raw files)
results/           tables and figures
manuscript/        the paper
```

## 7. Plan

| Step | Content | Estimate |
|---|---|---|
| 1 | Final variant sets: reuse the existing ClinVar code and the October count; add the gnomAD frequency flag | 2–3 days |
| 2 | Tool documentation and score retrieval (AVI, SpliceAI, Pangolin, CADD, phyloP) | 3–4 days |
| 3 | Part 1 analysis: tables and figures | 2–3 days |
| 4 | Part 2: obtain and harmonise the assay data, score, analyse | 1 week |
| 5 | Manuscript draft | 1 week |

Part 3 (prioritization in real genomes, including population representation) is decided after step 3. It can be a section of this paper or a separate paper.
