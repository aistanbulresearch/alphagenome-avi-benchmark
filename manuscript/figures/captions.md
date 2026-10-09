# Figure captions (draft)

**Figure 1. Discrimination on ClinVar labels added after the AlphaGenome training snapshot (Part 1).**
(a) ROC curves for intronic variants at least 3 bases from an exon in the main set (449 P/LP, 28,188 B/LB); false positive rate on a log scale. Pangolin was run on all P/LP and a region-stratified sample of 4,000 B/LB; its curve and AUROC use inverse-sampling weights. phyloP, 100-way vertebrate phyloP. (b) Paired AUROC differences, AVI minus each tool, on the variants both tools score, by region; intervals are 95% gene-cluster bootstrap intervals (2,000 replicates).

**Figure 2. Calibration as clinical evidence (Part 1).**
(a) Thresholds for each evidence level were derived on the calibration sample (labelled variants outside the new-label set) with the ClinGen method of Pejaver et al. (2022), prior 0.0441, and tested on the main set (non-canonical intronic and UTR variants): each level within its own score interval. Points are interval likelihood ratios with 95% gene-bootstrap intervals; vertical bars are the ClinGen targets (2.406, 5.790, 33.53 and 1,124; benign levels the reciprocals). Missing points: the level was not reached. Arrows: interval holding P/LP but no B/LB. (b) Share of P/LP receiving any pathogenic evidence and of B/LB receiving any benign evidence.

**Figure 3. Discrimination on experimental splicing assays (Part 2).**
AUROC difference, AVI minus (a) SpliceAI raw, (b) SpliceAI masked, (c) Pangolin, (d) CADD, per dataset (squares, 95% bootstrap intervals; variant bootstrap for single-gene datasets, gene bootstrap otherwise). Diamonds: random-effects pooled estimate over the main datasets (DerSimonian-Laird 95% interval); line through the diamond: Hartung-Knapp 95% interval. Grey: comparison datasets not pooled (MFASS; OpenSplice variants at the developers' benchmark positions). Pangolin was run on the four single-exon datasets.

**Figure 4. Part 2 by position relative to the exon.**
Pooled AUROC difference, AVI minus each tool, for exonic variants, intronic variants 3–50 bases from the exon, and intronic variants more than 50 bases from the exon (the Part 1 deep intronic cut). A dataset's class is included when it holds at least 10 splice-altering and 10 other variants; intronic >50 comprises COMPASS and OpenSplice only, so Pangolin (single-exon datasets) is absent there. Thick lines: DerSimonian-Laird 95% interval; thin lines: Hartung-Knapp 95% interval. Analysis added after the main results.

**Figure S1. SpliceAI2 (added comparator).**
(a) Part 1, AUROC difference AVI minus SpliceAI2 by region. (b) Part 2 per dataset and pooled over the five main datasets; OpenSplice, on which the SpliceAI2 developers evaluated their model, is shown separately and not pooled. Score: the SpliceAI2 summary score (precomputed v2.0).

**Figure S2. Missense CADD thresholds applied to intronic variants.**
The CADD score intervals recommended for missense variants (Pejaver et al. 2022, Table 2) evaluated as likelihood ratios on the intronic and UTR variants of the main set; bars are the ClinGen targets for each interval's evidence level; right, share of intronic P/LP falling in each interval. The intronic thresholds derived in this study are 17.6 (supporting), 19.2 (moderate) and 23.4 (strong). Analysis added after the main results.
