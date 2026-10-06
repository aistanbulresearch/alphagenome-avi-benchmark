# Data the AVI developers used, and how it is excluded

The new-label set must not contain variants that the AVI developers used for training, model selection or evaluation on the same task. This file lists those sources and how each is excluded. It is based on the AVI preprint: Methods, Data Split; AVI Evaluation Benchmarks; the GREGoR analysis; Data and Code Availability.

## Principle

A source's variants are excluded if either holds:

1. The source had a training or model-selection role.
2. The source was an evaluation benchmark measuring the same thing as our labels, that is pathogenicity or molecular/functional effect.

A source evaluated only on a different task (trait association) is not excluded.

## Sources

| Source | Role in the AVI preprint | Measures | Exclusion |
|---|---|---|---|
| Saturation genome editing (SGE) | Model selection on the validation chromosomes (four studies); test (ten studies) | Functional effect | Positions of every scored variant, or the gene region (below) |
| TraitGym Mendelian | Test | Pathogenicity (OMIM) | Positions of every variant |
| TraitGym Complex | Model selection on the validation chromosomes; test | Fine-mapped trait association | Positions of every variant |
| ncVarDB | Test | Pathogenicity | Positions of every variant |
| GREGoR | Test (retrospective ranking of 112 P/LP variants) | Pathogenicity | Submitter rule (below); the variant list is not published |
| Complex-trait benchmark | Test only | Trait association | Not excluded; the variant list is not published |

The preprint gives the chromosome groups as follows. Validation is SGE and TraitGym Complex variants on chromosomes 2, 5, 11, 14, 17, 20, 22 and X. Test is positions not used in training or validation.

## Saturation genome editing

Every GRCh38 position of every variant scored in the listed files is excluded, whatever its functional class, using the full span for indels.

- **MaveDB studies:** scores come from `https://api.mavedb.org/api/v1/score-sets/<urn>/scores`. Positions come from MaveDB's mapped variants for the same score set, covering every version of a listed score set.
- **Gene-region fallback:** if a file cannot be obtained, or a scored variant has no GRCh38 position, the study's whole gene is excluded instead.
- **What the gene region is:** the GRCh38 locus of the GeneID, padded by 1,000 bases on each side. It is taken from the `gene` feature of NCBI RefSeq annotation release GCF_000001405.40-RS_2025_08 (the release MANE v1.5 names), using one span per chromosome.
- **Studies without a MaveDB score set** use the gene region: De Jonghe 2026 (RNU4-2), Lee 2025 (ATM) and Dace 2025 (BRCA1).

| Study | Gene | Files |
|---|---|---|
| Findlay et al. 2018, Nature 562:217 | BRCA1 | MaveDB `urn:mavedb:00000097-0` and exon score sets `00000097-a` to `-z` |
| Huang et al. 2025, Nature 638:528 | BRCA2 | `urn:mavedb:00001225-a`, `urn:mavedb:00001279-a` |
| Waters et al. 2024, Nat Genet 56:1434 | BAP1 | `urn:mavedb:00000662-0` and `00000662-a` to `-r` |
| Olvera-León et al. 2024, Cell 187:5719 | RAD51C | `urn:mavedb:00000673-0` and `00000673-a` to `-i` |
| Buckley et al. 2024, Nat Genet 56:1446 | VHL | `urn:mavedb:00000675-a` |
| Radford et al. 2023, Nat Commun 14:7702 | DDX3X | `urn:mavedb:00000658-0` and `00000658-a` to `-q`; `urn:mavedb:00001228-a` |
| De Jonghe et al. 2026, Nature 654:429 | RNU4-2 | gene region |
| Lee et al. 2025, Cell | ATM | gene region |
| Dace et al. 2025, preprint | BRCA1 | gene region |
| Woo et al. 2025, preprint | BARD1 | `urn:mavedb:00001250-a` |
| Tejura et al. 2026, bioRxiv | PALB2 | `urn:mavedb:00001259-a` |
| Tejura et al. 2026, bioRxiv | XRCC2 | `urn:mavedb:00001264-a` |
| Tejura et al. 2026, bioRxiv | SFPQ | `urn:mavedb:00001265-a` |
| Tejura et al. 2026, bioRxiv | CTCF | `urn:mavedb:00001262-a` |

In the October 2026 build every study fell back to its gene region, which is the conservative choice. The genes excluded this way are reported in the results.

## TraitGym and ncVarDB

| Source | Version | Files |
|---|---|---|
| TraitGym (Mendelian and Complex) | Hugging Face `songlab/TraitGym`, revision `1fde19555fe8c0a55b1382bdf4c6f7082209f566` (2025-03-25) | `mendelian_traits_all`, `mendelian_traits_matched_9`, `mendelian_traits_v21_matched_9`, `mendelian_traits_v22_matched_9`, `complex_traits_all`, `complex_traits_matched_9`, `complex_traits_v22_matched_9` (each `test.parquet`) |
| ncVarDB | GitHub `Gardner-BinfLab/ncVarDB`, commit `42828402cd3fea565308faaefa187189faf2ec79` | `data/ncVar_pathogenic.tsv`, `data/ncVar_benign.tsv`, `data/ncVar_pathogenic.vcf`, `data/ncVar_benign.vcf` |

The positions of every variant are excluded, positives and controls alike.

## GREGoR: submitter rule

The preprint takes its 112 variants from the Broad Institute Center for Mendelian Genomics subset of GREGoR. That variant list is not published.

- **The rule:** every ClinVar variant with a submission from ClinVar organisation 506627 ("Broad Center for Mendelian Genomics; Broad Institute of MIT and Harvard") is excluded. Both spellings of the name found in ClinVar are matched.
- **What it misses:** a variant among the 112 that only another organisation submitted remains in the set.

## Complex-trait benchmark: not excluded

The developers evaluated AVI on rare variants associated with complex traits in UK Biobank and All of Us. That is a different task from pathogenicity, and the variant list is not published.
