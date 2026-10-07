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
  - Variants at positions in benchmark datasets the developers used for training, selection or evaluation (saturation genome editing screens, TraitGym, ncVarDB; see `data/development_sources.md`).
  - Submissions from the organisation that supplied the developers' GREGoR cases.
  - Variants at the position of an older strict ClinVar record of at most 50 bases whose identity cannot be linked to any current record (it may be the same variant under an older identifier).
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

### 4.4 Analysis plan, fixed before any score is retrieved

**Main analysis set.** New labels; SNVs; at least one review star; the variant's HGVS name is on its gene's MANE Select transcript; untouched and validation chromosomes; rare in gnomAD v4.1.

- **Rare:** the joint filtering allele frequency (`fafmax_faf95_max_joint`) is below 0.001, the threshold AVI's training used. If the variant is absent from gnomAD, or the field is empty, the joint allele frequency below 0.001 is used, and absence counts as rare.
- **Why rare is the main set:** a variant that is common in the population is classified by its frequency, not by a predictor, and frequency is what AVI was trained on.

**Supporting analyses**, each reported in full next to the main one:

- all frequencies;
- rare variants present in gnomAD only (absence can reflect low coverage);
- training chromosomes;
- at least two stars;
- the stricter exclusion of variants inside longer unknown older records (column `inside_long_unknown_record`), so that the result does not depend on that data-cleaning choice;
- by whether any ClinVar submission of the variant mentions a splicing predictor (`data/derived/variants/submission_flags.tsv.gz`; presence of fixed terms only, no values), because labels that rest on SpliceAI favour SpliceAI;
- deep and near intronic SNVs scored with SpliceAI run locally with a ±4,999-base window (pseudoexon-creating variants can lie beyond the 50-base window of the precomputed scores), if the local run is feasible;
- the calibration sample (below), which represents the labelled variants outside the new-label set;
- small indels, for tools other than AVI (AVI publishes SNV scores only).

**Regions.**

- **Main comparison region:** intronic variants at least 3 bases from the exon, which pools the donor region, acceptor region, near intronic and deep intronic.
- **Reported separately:** each of those four regions, canonical sites (±1, ±2), UTR (5' and 3' pooled, and each alone) and synonymous.

**Scores.**

| Tool | Score used | Source |
|---|---|---|
| AVI | raw score of the `AVI_SCORE` scorer, with the calibrated score the API returns | AlphaGenome Atlas API (same precomputed scores as the AVI SNV download) |
| SpliceAI | maximum of the four delta scores, raw (unmasked), distance 50; the maximum over gene entries. Masked scores as a secondary form | Ensembl MANE v1.4 precomputed files |
| CADD v1.7 | PHRED score | whole-genome SNV file |
| Pangolin | maximum absolute splice score change, masked, distance 50 | local run |
| phyloP | 100-way vertebrate, per position | UCSC bigWig |
| AVI's conservation inputs | phastCons 470-way; Cactus 241-way phyloP | UCSC bigWig |
| Distance baseline | distance to the nearest exon boundary (from the HGVS offset), smaller = more likely pathogenic | variant list |

- **Missing scores:** a variant a tool cannot score is counted in a coverage table. Comparisons between two tools use the variants both tools score.
- **One retrieval:** the scores of every set (new labels and calibration sample) are retrieved in one pass, so no set is shaped after another set's scores are seen.
- **AVI source check:** before the main retrieval, the API's raw values are compared with the downloadable AVI file on a fixed sample of positions outside the study sets.
- **AVI's own inputs:** AVI against phastCons 470-way and Cactus 241-way phyloP is reported as "what AVI adds beyond its conservation inputs".

**Discrimination.**

- **Metrics:** AUROC and AUPRC per tool, set and region. The P/LP fraction, which is AUPRC's baseline, is reported with each AUPRC.
- **Main comparison:** the paired AUROC difference, AVI minus SpliceAI, in the main region of the main set. AVI is also compared with CADD in every region, and with SpliceAI in every intronic region.
- **Intervals:** 95% intervals from 2,000 bootstrap replicates that resample genes (seed 20261007). Every difference is reported with its interval; no difference is selected for reporting by its size.

**Clinical evidence strength.**

- **Method:** the ClinGen calibration of Pejaver et al. (2022), using its local posterior probability and a prior of 0.0441. Score thresholds are derived for supporting, moderate and strong pathogenic evidence, and for the corresponding benign levels. The posterior targets are those of Pejaver et al. 2022, Table 1.
- **Where thresholds come from:** they are derived on a calibration sample of the labelled variants outside the new-label set. It contains non-canonical intronic and UTR variants with at least one star, then restricted to rare ones. Every P/LP variant is included, plus 40,000 B/LB variants drawn at random, stratified by region and chromosome group (seed 20261007; list in `data/derived/variants/calibration_sample.tsv`, committed before any score). The calibration weights the classes to the prior, so the B/LB sampling fraction does not change the posterior.
- **Not calibrated:** canonical ±1/±2 sites, which are evidence for a loss-of-function criterion rather than a predictor criterion, and synonymous variants.
- **How they are tested:** on the main set, as the positive and negative likelihood ratios at each threshold, with gene-bootstrap 95% intervals.
- **Published thresholds:** SpliceAI's ClinGen thresholds (Walker et al. 2023: at least 0.2 for PP3, at most 0.1 for BP4) are tested on the main set. Walker et al. defined them on the raw maximum delta score with a ±4,999-base window; our main SpliceAI score is raw with a 50-base window, and the difference is stated.

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
