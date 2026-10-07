# Part 2 sources: experimental splicing measurements

Built by `src/part2_data.py` (`download`, then `build`) into `data/derived/part2/assay_variants.tsv.gz`.
Every download is recorded in `data/downloads.json` (URL, date, bytes, SHA-256). Raw files are not in the
repository; they are read from the raw folder (`AVI_RAW_DIR`) under `part2/`.

No tool score was retrieved, computed or read while building this table. Source tables that also carry
prediction columns are read with an explicit list of assay columns, so prediction columns are never loaded;
supplementary files that contain only predictions were not downloaded.

## Common conventions

- **Coordinates.** GRCh38, chromosome without `chr`, 1-based position, forward-strand alleles. Indels are
  left-aligned in VCF form (one anchor base). hg19 positions are converted with `src/liftover.py` and the UCSC
  chain `data/derived/exposure/hg19ToHg38.over.chain.gz`.
- **Reference check.** Every reference allele is compared with the GRCh38 FASTA (`hg38.fa`, read through its
  `.fai`). Rows that do not match are dropped and counted (table below).
- **Labels.** `splice_altering` is the authors' own call, or the authors' own threshold applied to the
  authors' own effect size where the published table has no call column. `intermediate` is used only where
  the authors define it. `effect_size` keeps the continuous measurement named in `effect_measure`.
- **Use.** `main` = selection-free library variants; `comparison` = MFASS (developers of compared tools
  report results on it); `not_used` = subsets excluded by design, duplicates, or uninformative measurements
  (reason in `subset`).
- **Duplicates.** If the same allele occurs more than once within a dataset among usable rows, the first row in
  `variant_id_in_source` order is kept and the others get `use = not_used` and `;duplicate` in `subset`. The
  rule does not look at labels or effect sizes.
- **Variant classes.** `variant_type` is SNV, deletion, insertion or delins. Indels stay in the table with
  their own type; AVI publishes SNV scores only.

## Counts

Rows = rows read from the source table (variants). Usable = `use` main or comparison. Labels are for usable
SNVs (yes / no / intermediate).

| Dataset (`dataset`) | Source rows | Ref mismatch dropped | In table | Usable SNVs | SNV yes / no / int. | Usable indels (yes / no) | Use |
|---|---|---|---|---|---|---|---|
| POU1F1 exon 2 (`kitzman_pou1f1_ex2`) | 1,070 | 0 | 1,070 | 1,070 | 96 / 938 / 36 | 0 | main |
| WT1 exon 9 (`kitzman_wt1_ex9`) | 518 | 0 | 518 | 518 | 57 / 445 / 16 | 0 | main |
| MST1R exon 11 (`braun2018_mst1r_ex11`) | 2,029 | 3 | 2,026 | 1,850 | 876 / 974 / 0 | 173 (122 / 51) | main |
| CD19 exons 1-3 (`cortes_lopez2022_cd19`) | 4,255 | 4 | 4,251 | 3,596 | 178 / 3,418 / 0 | 655 (11 / 644) | main |
| COMPASS (`compass`) | 87,546 | 0 | 35,712 | 32,246 | 7,138 / 25,108 / 0 | 0 | main (random SNV, ExAC) |
| OpenSplice (`opensplice`) | 882,382 | 0 | 713,622 | 310,369 | 171,983 / 138,386 / 0 | 286,952 (190,210 / 96,742) | main |
| MFASS (`mfass`) | 32,669 | 91 | 27,642 | 27,642 | 1,048 / 26,594 / 0 | 0 | comparison |

Details per dataset are in the sections below. The build prints the full count report (`--summary`).

---

## 1. POU1F1 exon 2 saturation splicing assay (Kitzman laboratory)

- **Citation.** Gergics P, Smith C, Bando H, et al. High-throughput splicing assays identify missense and
  silent splice-disruptive POU1F1 variants underlying pituitary hormone deficiency. Am J Hum Genet
  2021;108(8):1526-1539. doi:10.1016/j.ajhg.2021.06.013 (PMC8387473).
- **Data.** Zenodo 10.5281/zenodo.7803918, file `POU1F1_SpliceEffects.xls` (Table S1 of the paper;
  `FormatManual.xlsx` describes the record). Licence **CC-BY-4.0** (Zenodo record).
- **Assay.** Pooled exon-trap minigene (pSPL3) with POU1F1 exon 2 (beta) and flanking intron; RNA-seq of
  isoforms (alpha, beta, skip, other); replicate transfections (calls require at least half the replicates).
- **Context.** Construct (minigene), not the endogenous locus.
- **Variant choice.** Every possible SNV across exon 2 beta (150 bp) and 210 bp of flanking intron; 1,070 of
  1,080 designed SNVs were measured. Not chosen by any predictor.
- **Threshold (authors).** Section "Fold-change and significance testing": SDV = Z > 4.16 (Bonferroni,
  p = 0.05) and "at least a 3-fold change from the average null distribution PSI" for the beta, skip or other
  isoform in at least half the replicates; intermediate = Z met with a 2- to 3-fold change, or SDV criteria met
  overall but not individually in at least 7 replicates. The table's `Splice Disruptive` and `Intermediate`
  columns are used as published.
- **Effect size.** Maximum of the beta, skip and other isoform fold changes (the authors' fold-change metric).
- **Cell type.** COS-7 (African green monkey kidney cells).
- **Mapping.** The table gives GRCh38 positions; alleles are on the transcript (minus) strand and are
  complemented (all 1,070 then match GRCh38; 0 of 1,070 matched as given).
- **Exposure.** No training or evaluation use by the developers of AVI, AlphaGenome, SpliceAI, Pangolin or
  CADD was found (dataset inventory of 2026-10-06). A third-party benchmark (Smith & Kitzman 2023, Genome Biol,
  doi:10.1186/s13059-023-03144-z) evaluated several splicing predictors, including SpliceAI and Pangolin, on it.

## 2. WT1 exon 9 saturation splicing assay (Kitzman laboratory)

- **Citation.** Smith C, Burugula BB, Dunn I, Aradhya S, Kitzman JO, Yee JL. High-throughput splicing assays
  identify known and novel WT1 exon 9 variants in nephrotic syndrome. Kidney Int Rep 2023;8(10):2117-2125.
  doi:10.1016/j.ekir.2023.07.033 (PMC10577367).
- **Data.** Zenodo 10.5281/zenodo.7803918, file `WT1_SpliceEfects.xlsx` (Supplementary Table 2). Licence
  **CC-BY-4.0** (Zenodo record; the article itself is CC BY-NC-ND).
- **Assay.** pSPL3 minigene with WT1 exon 9 and 414 bp of flanking intron; isoforms KTS+, KTS-, SKIP, OTHER.
- **Context.** Construct.
- **Variant choice.** All possible SNVs in the exon and 40 bp into each intron (518 of 519 measured). Not chosen
  by any predictor.
- **Threshold (authors).** Methods, "RNA-seq Processing and Splice Disruption Calling": SDV = Bonferroni
  P < 0.05 and either SKIP or OTHER usage at least 20% above the null, or "an isoform log-ratio (calculated as
  log2[KTS+/KTS−]) of ≥1.5 or ≤ −1"; intermediate = same test with 10% above the null or a log-ratio
  >= 1 or <= -0.5. Published `Splice Disruptive` and `Intermediate` columns are used.
- **Effect size.** Signed log2(KTS+/KTS-) ratio. SKIP and OTHER usage, which drive a minority of the SDVs
  (14 exon skipping and 2 cryptic acceptor of 57), are not in this single number.
- **Cell type.** HEK293T (8 replicates; 2 COS-7 replicates were also run).
- **Mapping.** HGVS c. alleles on the transcript (minus) strand, complemented; hg19 positions lifted to GRCh38
  (518 of 518 lifted and matched).
- **Exposure.** As for POU1F1 (third-party benchmark only).

## 3. MST1R (RON) exon 11 random-mutagenesis minigene

- **Citation.** Braun S, Enculescu M, Setty ST, et al. Decoding a cancer-relevant splicing decision in the RON
  proto-oncogene using high-throughput mutagenesis. Nat Commun 2018;9:3315. doi:10.1038/s41467-018-05748-7
  (PMC6098099).
- **Data.** Supplementary Data 3 (`41467_2018_5748_MOESM6_ESM.xlsx`) and the description of supplementary files
  (`MOESM3`). Licence **CC BY 4.0** (article and supplementary material).
- **Assay.** Minigene with exon 11, its complete flanking introns and constitutive exons 10 and 12, mutated by
  error-prone PCR (target 3.5 mutations per minigene; 5,200 mutant and 591 wild-type barcoded minigenes);
  targeted RNA-seq of six isoform classes, 3 replicates.
- **Context.** Construct.
- **Variant choice.** Random mutagenesis; not chosen by any predictor.
- **Modelling caveat.** Single-mutation effects are not measured directly: they are inferred by a linear model of
  log fold changes over multi-mutant minigenes (Supplementary Note 2). 131 mutations that the model could not
  infer are given only as median-based estimates.
- **Threshold (authors).** Methods, "Significant mutation effects and synergistic interactions": a mutation is
  significant for an isoform with "≥5% change in isoform frequency compared to the mean of the 591 minigene
  variants" and FDR < 5% (Stouffer-combined replicate z-scores, Benjamini-Hochberg). The published table has no
  call column, so the rule is applied here to the published replicate estimates and wild-type means and SDs
  (z per replicate, Stouffer, two-sided p, BH within each isoform, change = mean over replicates). This
  reproduces the paper's 778 splicing-effective point mutations in HEK293T exactly (775 remain after 3 reference
  mismatches). Label = significant for any of the six isoforms.
- **Median-based estimates** (`subset = median_estimate`, 106 SNVs, 25 indels): no test can be applied; labelled
  by |median change| >= 5 points for any isoform. Kept as `main`, separable by `subset`.
- **Effect size.** Maximum |change in isoform frequency| over the six isoforms (percentage points).
- **Cell type.** HEK293T (the MCF7 and knockdown measurements are not used).
- **Mapping.** hg19, minus strand; hg19 position = 49,933,838 - minigene position (one constant for every SNV).
  For indels the anchor base is taken from the authors' genomic position, because the number in the indel
  notation is the anchor in one sheet and the deleted base in the other. Lifted to GRCh38. 3 SNVs at minigene
  position 298 (G298A/C/T) do not match GRCh38 and are dropped.
- **Exposure.** None found for the developers of the compared tools; third-party benchmark (Smith & Kitzman
  2023).

## 4. CD19 exons 1-3 random-mutagenesis minigene

- **Citation.** Cortés-López M, Schulz L, Enculescu M, et al. High-throughput mutagenesis identifies mutations
  and RNA-binding proteins controlling CD19 splicing and CART-19 therapy resistance. Nat Commun 2022;13:5570.
  doi:10.1038/s41467-022-31818-y (PMC9500061).
- **Data.** Supplementary Data 3 (`41467_2022_31818_MOESM6_ESM.xlsx`, sheets "Mutation effects" and "19 validated
  mutations"), Supplementary Data 2 (`MOESM5`), description file (`MOESM3`). Licence **CC BY 4.0**. Raw reads:
  GEO (series listed in the paper's Data availability statement). Supplementary Data 5 (predictions only) was not
  downloaded.
- **Assay.** Minigene amplified from genomic DNA (CD19 exons 1-3 with introns), error-prone PCR (9.7 mutations
  per minigene on average), 9,321 minigene variants; targeted RNA-seq, 2 replicates; single-mutation effects on
  six isoform classes inferred by a penalised multinomial (softmax) regression.
- **Context.** Construct.
- **Variant choice.** Random mutagenesis; not chosen by any predictor.
- **Threshold (authors).** Results, "Therapy-impacting isoforms accumulate in response to numerous point
  mutations", and Methods: a mutation affects an isoform if, in both replicates and in the same direction, its
  predicted frequency lies "beyond the 2.5 and 97.5% quantiles of the WT minigene distribution". The published
  `splicing-effective mutation` column is used (Yes for any isoform). 193 mutations are called in the paper; 189
  remain after 4 reference mismatches.
- **Effect size.** Maximum |average change in isoform frequency| over the six isoform classes (fraction).
- **Cell type.** NALM-6 (B-cell precursor leukaemia).
- **Mapping.** GRCh38, plus strand: position = 28,931,870 + minigene position, derived from the authors' own
  hg38 coordinates of validated mutations (T465A, C645T; one constant). Indel notation includes the anchor
  base. 4 mutations at minigene positions 747-748 (GT747G, T748A/C/G) do not match GRCh38 (the minigene base
  differs from the reference) and are dropped.
- **Exposure.** None found for the developers of the compared tools. The authors ran SpliceAI on the CD19 gene
  (third-party use).

## 5. COMPASS (Seelig laboratory)

- **Citation.** Koplik SE, Yu AM, Shelby MR, et al. Massively parallel assay of human splice variants reveals
  cis-regulatory drivers of disease-associated variants (title of the PMC record). bioRxiv 2025.
  doi:10.1101/2025.10.12.681955 (PMC12633039; GEO series title "Deciphering the cis-regulatory code of
  alternative splicing in human cell lines").
- **Data.** GEO GSE307247: `GSE307247_Processed_PSIs_All_Cells.csv.gz`, `GSE307247_column_descriptions_metadata.tsv.gz`,
  series SOFT file. Licence: GEO gives no licence statement for the processed file (NCBI places no restriction on
  the use of GEO data but does not grant rights on submitter content); preprint **CC BY-NC 4.0**; code repository
  github.com/skoplik/ESL_MPRA_2025 **MIT** (design table `st_final_with_snp_and_coords_05_30_25.csv.gz` at commit
  67791aa was consulted for the source labels; it does not separate ClinVar from random variants either).
- **Assay.** Exon-skipping MPRA: 161-bp variable region (exon of at most 90 nt, 20 nt of downstream intron and
  upstream intron) between constant SMN2 intron 6/7 sequences in a Citrine reporter; PSI from pooled replicates.
- **Context.** Construct.
- **Variant choice.** 2,096 short exons from 1,733 genes. Subsets (authors' labels): references; ExAC variants;
  ClinVar variants; Geuvadis variants; random SNVs and random double variants "randomly targeted" to introns and
  exons containing at least one ClinVar variant (region choice by ClinVar, variant choice random).
  - **main:** random SNVs (`random_snv`, 28,477) and ExAC SNVs (`exac`, 3,769).
  - **not used:** ClinVar subset (`clinvar`, 1,885), Geuvadis subset (`geuvadis`, 47); doubles (48,290 sequences)
    and references (1,323) are not in the table. The published source label `clinvar_single` mixes the ClinVar
    variants and the random SNVs; an allele of that label present in ClinVar 2025-06-15 (the developers'
    snapshot, already in the raw folder) is assigned to the ClinVar subset, the rest to the random SNVs.
- **Threshold (authors).** Results, "Quantifying and predicting variant impact on splicing": "we considered any
  variant with ∣Δlogit(PSI)∣ ≥ 1 to have a significant effect on splicing".
- **Effect size.** Mean over cell lines of the pooled-replicate delta logit(PSI), as in the authors' global
  analyses (Methods, "Bioinformatic Workflow"); the cell lines measured are listed in `effect_measure`. The
  processed table already applies the authors' coverage filter (at least 10 reads per replicate).
- **Uninformative measurements.** The authors exclude variants whose reference PSI is exactly 0 or 1 when they
  evaluate predictors (Results, same section). Main-set variants whose reference PSI is 0 or 1 in every cell line
  where they were measured are `not_used` (`;reference_psi_0_or_1`, 1,426). 2,221 single variants have no
  measurement in any cell line and are not in the table (2,122 of them from the ExAC subset).
- **Cell types.** HEK293, HeLa, K562, MCF7, HMC3.
- **Mapping.** `variant_hg38` (GRCh38, forward strand) is used as published; all 35,712 match GRCh38. 108 alleles
  that occur in two sequences are flagged as duplicates.
- **Exposure.** No developer use found for AVI (the AVI preprint, 2026-09-16, does not mention COMPASS) or the
  other compared tools. The COMPASS authors compared their measurements with several predictors, including
  SpliceAI and Pangolin, and their code repository contains AlphaGenome prediction scripts (third-party use). ExAC variants may coincide with AVI's gnomAD-derived
  training variants; positional overlap is to be checked when scores are joined.

## 6. OpenSplice (Lehner and Valcárcel laboratories)

- **Citation.** Quarantani G, Clarke J, Thompson M, Sang F, Valcárcel J, Lehner B. OpenSplice: the impact of half
  a million mutations on the alternative splicing of 600 human exons. bioRxiv 2026 (posted 2026-05-23).
  doi:10.64898/2026.05.22.727141.
- **Data.** Publicly downloadable. Per-variant PSI = Supplementary Table 4 of the preprint (bioRxiv supplementary
  file `media-6.zip`, `Supplementary_Table4.tsv`). Licence: the preprint and its supplementary material are
  **CC BY-NC 4.0** (bioRxiv). Related: Figshare 10.6084/m9.figshare.32337414 (CC BY 4.0; barcode-level counts
  and predictor files, not downloaded), ENA PRJEB111846 (reads), code github.com/lehner-lab/OpenSplice (MIT).
  Supplementary Table 12 (predictions) was not downloaded.
- **Assay.** Site-saturation mutagenesis of 608 human exons with proximal flanking intron, each inserted into a
  FAS exons 5-7 minigene in place of FAS exon 6; RNA-seq of inclusion and skipping, 3 replicates; PSI from a
  hierarchical error model.
- **Context.** Construct (non-native FAS flanks).
- **Variant choice.** All SNVs and all 1, 3, 6 and 21 nt deletions in each tested exon and its proximal introns
  (not chosen by any predictor). Exons were chosen from a screen of 6,000 wild-type exons, prioritising
  intermediate inclusion, alternative splicing and the number of pathogenic ClinVar variants (ClinVar
  2023-07-22) per exon, i.e. region choice partly by clinical annotation, variant choice exhaustive.
- **Threshold (authors).** Methods, "PSI calculation / Variant effect testing": "A variant was called a
  significant splicing-altering when the adjusted p-value < 0.1" (z-test on the logit scale, BH within each
  exon). The published `significant` column is used.
- **Effect size.** `delta_logit` (mutant minus wild-type logit PSI).
- **Cell type.** HEK293T.
- **Mapping.** CHROM/POS/REF/ALT are GRCh38; deletions are given in VCF form with a "start_enddel" position and
  are placed on their anchor base and left-aligned. All 713,622 measured rows match GRCh38. 112,332 wild-type
  rows and 56,428 rows without a measurement are not in the table. 116,301 duplicates (116,265 deletions that
  become the same allele after left alignment, and 36 SNVs measured in two exon constructs, 14 of them with
  different calls) are `not_used`. Usable counts (310,369 SNVs; 286,952 deletions) agree with the paper's
  310,405 substitutions and 287,947 deletions.
- **Exposure.** The OpenSplice authors benchmarked SpliceAI, Pangolin, AlphaGenome and another predictor as third
  parties. No developer use found; the dataset was posted before the AVI preprint. Several tested exons lie in
  genes whose saturation genome editing screens were AVI development material (ATM, BAP1, BARD1, BRCA1, BRCA2,
  PALB2, RAD51C; 26,080 usable rows).

## 7. MFASS (comparison only)

- **Citation.** Chong R, Insigne KD, Yao D, et al. A multiplexed assay for exon recognition reveals that an
  unappreciated fraction of rare genetic variants cause large-effect splicing disruptions. Mol Cell
  2019;73(1):183-194.e8. doi:10.1016/j.molcel.2018.10.037 (PMC6599603).
- **Data.** github.com/KosuriLab/MFASS, `processed_data/snv/snv_data_clean.txt` at commit 9a8e4f2. Licence:
  **none stated** in the repository.
- **Assay.** Sort-seq minigene assay (MFASS): 170-bp exon-plus-intron library sequences in a reporter integrated
  as a single copy at AAVS1; exon inclusion index from sorted bins (library version 2 used in the paper).
- **Context.** Construct.
- **Variant choice.** All ExAC PASS SNVs within nucleotides 11-160 of each library sequence, 2,198 exons
  (Ensembl 73, hg19 design). Selection-free with respect to predictors, but used as a comparison set only
  (below).
- **Threshold (authors).** Results, "Evaluating MFASS based on Known Splicing Regulatory Elements": "we defined
  large-effect variants as Δinclusion index ≤ −0.5", for a wild-type exon with inclusion index >= 0.5 (SDV). The
  published `strong_lof` column is used; 1,239 mutant rows without a call are excluded (27,733 called, 1,050 SDVs,
  as in the paper).
- **Effect size.** `v2_dpsi` (delta inclusion index, library v2).
- **Cell type.** HEK293T.
- **Mapping.** The authors' GRCh38 positions (`snp_position_hg38_1based`) and forward-strand alleles. 91 rows whose
  reference allele matches GRCh38 only as the reverse complement are dropped (2 of them SDVs).
- **Exposure.** AlphaGenome (Avsec et al., Nature 2026) and Pangolin (Zeng & Li, Genome Biol 2022) report
  evaluations on MFASS; CADD-Splice used MFASS when choosing its splicing features (dataset inventory; not
  re-verified). Hence `use = comparison`.

## Overlaps noted at build time (usable SNVs, identities only)

- Between datasets: COMPASS x OpenSplice 1,806; COMPASS x MFASS 118; MFASS x OpenSplice 129.
- With the part 1 ClinVar list (`data/derived/variants/clinvar_variants.tsv.gz`), all labelled / new labels:
  POU1F1 20 / 0; WT1 64 / 8; MST1R 6 / 1; CD19 48 / 2; COMPASS 229 / 114; OpenSplice 12,845 / 971;
  MFASS 2,293 / 178. STUDY_DESIGN.md section 5.1 counts a variant present in both parts in one part only; that
  assignment is not made in this table.
