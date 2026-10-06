# AlphaGenome AVI benchmark

An independent evaluation of the AlphaGenome Variant Impact (AVI) score for non-coding and splice-affecting variants. AVI is compared with established tools on two kinds of ground truth:

1. ClinVar pathogenic and benign classifications added after the ClinVar release used by the AVI developers.
2. Experimental splicing assays whose variants were not selected by any tool's prediction.

The study design is in [STUDY_DESIGN.md](STUDY_DESIGN.md).

## Repository layout

| Path | Content |
|---|---|
| `STUDY_DESIGN.md` | question, variant sets, tools, metrics, plan |
| `src/` | code |
| `tests/` | tests of the code |
| `data/` | list of downloaded sources (URL, date, checksum) and derived variant lists; raw downloads are not stored here |
| `results/` | tables and figures |
| `manuscript/` | the paper |

## Reproducing

Instructions are added as each step is completed: the download list, the commands that build the variant sets, retrieve the scores and compute every table and figure.

## Licence

Code: MIT (see `LICENSE`). Third-party data and tool scores keep their own terms; files that may not be redistributed are not included, and the download list says how to obtain them.
