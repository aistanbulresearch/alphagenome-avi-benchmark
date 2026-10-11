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
| `data/` | list of downloaded sources (URL, date, checksum), source records of the score files, and derived variant lists without scores; raw downloads are not stored here |
| `results/` | tables and figures |
| `manuscript/` | figures and legends of the paper; `manuscript/tools/` builds the numbers sheet and the supplementary tables |

## Reproducing

Every number, table and figure of the manuscript is computed by the scripts in `src/` from public files. Commands are run from the repository root. Raw downloads are not stored here: `data/downloads.json` lists each file with its URL, date and SHA-256, and `data/derived/scores/*_source.json` records where each score file came from: the source file or API, the version, the retrieval or run time and the variant counts, with the checksum of the remote file's index where one was read, and for the local runs the package, environment, genome checksum and settings. Per-variant tool scores are not distributed (the SpliceAI2 terms forbid it; for the other tools only the sources are given). Variant lists without scores are in `data/derived/variants/`, including the Part 1 main set (`main_set.tsv`).

### Environments, folders and keys

- Analysis: Python 3.14 with `pip install -r requirements.txt` (the AlphaGenome client is pinned to 0.9.0).
- SpliceAI and Pangolin run locally in a second environment: Python 3.11, SpliceAI 1.3.1 with tensorflow-cpu 2.21, and Pangolin 1.0.2 from github.com/tkzeng/Pangolin at commit 5cf94b8 with torch 2.14.1 (CPU), installed without its VCF dependencies (versions in `data/derived/scores/spliceai_d4999_source.json` and `pangolin_source.json`).
- Raw files go to `data/raw/`, or to the folder named by the environment variable `AVI_RAW_DIR`.
- `ALPHAGENOME_API_KEY` (AlphaGenome Atlas API, for AVI) and `HF_TOKEN` (the SpliceAI2 dataset on Hugging Face, after accepting its terms) are read from the environment and never written to a file.

### 1. Downloads

For each entry of `data/downloads.json`, `python src/download_public.py --name <path> <url>`, where `<path>` is the entry's `file` without the leading `data/raw/`; check the SHA-256 against the recorded one (the script appends a new record to the list). The MaveDB files are fetched with `python src/mavedb_sources.py fetch` and the Part 2 sources with `python src/part2_data.py download`. Decompress `hg38.fa.gz` to `hg38.fa` for the local SpliceAI, Pangolin and Part 2 runs.

### 2. Exposure files (benchmark and training sources of the AVI developers)

```
python src/build_exposure_files.py submitters
python src/build_exposure_files.py chain
python src/build_mane_table.py
python src/build_replacements.py
python src/build_development_positions.py
```

`data/development_sources.md` describes the sources; the outputs other than `development_positions.tsv` and the chain file are kept in `data/derived/exposure/`.

### 3. Part 1 variant list (ClinVar)

```
python src/submission_flags.py --submission-summary data/raw/submission_summary.txt.gz
python src/variant_sets.py --release data/raw/variant_summary_2026-10.txt.gz --e-dev data/raw/clinvar_20250615.vcf.gz --t0 data/raw/variant_summary_2025-07.txt.gz --e-hist data/raw/variant_summary_2024-07.txt.gz --submission-summary data/raw/submission_summary.txt.gz
python src/gnomad_frequencies.py
python src/calibration_sample.py
python src/gnomad_frequencies.py --rows calibration
```

### 4. Part 1 scores

```
python src/avi_atlas.py --chunk 500 --workers 20
python src/tool_scores.py --tool spliceai --threads 16
python src/tool_scores.py --tool spliceai_masked --threads 16
python src/tool_scores.py --tool cadd --threads 16
python src/conservation_scores.py --threads 16
```

In the second environment:

```
python src/spliceai_local.py --fasta data/raw/hg38.fa --batch 32
python src/pangolin_local.py --fasta data/raw/hg38.fa --gtf data/raw/gencode.v38.annotation.gtf.gz --batch 16
```

The Part 1 AVI scores come from `src/avi_atlas.py` (Atlas API, client 0.9.0; `avi_source.json`). The `--tool avi` path of `src/tool_scores.py`, which reads a local copy of the downloadable AVI file, was not used for the reported results. The SpliceAI2 scores of both parts are retrieved in step 5, because `src/spliceai2_scores.py` also reads the Part 2 alleles.

### 5. Part 2 assay data and scores

```
python src/part2_data.py build
python src/part2_scores.py --tool avi --threads 16
python src/part2_scores.py --tool spliceai --threads 16
python src/part2_scores.py --tool spliceai_masked --threads 16
python src/part2_scores.py --tool cadd --threads 16
python src/part2_scores.py --tool conservation --threads 16
python src/spliceai2_scores.py
```

In the second environment: `python src/pangolin_local.py --fasta data/raw/hg38.fa --gtf data/raw/gencode.v38.annotation.gtf.gz --batch 16 --part2`.

### 6. Analyses, figures and manuscript tables

```
python src/part1_flow.py
python src/part1_analysis.py
python src/part1_supplement.py
python src/part2_analysis.py
python src/part2_by_position.py
python src/part2_positions.py
python src/part2_compass_subsets.py
python src/part2_compass_clinvar_check.py
python src/export_main_set.py
python src/figures.py
python manuscript/tools/numbers_sheet.py > NUMBERS.md
python manuscript/tools/make_supplementary_tables.py supplementary_tables.xlsx
```

Results are written to `results/part1/` and `results/part2/`, figures to `manuscript/figures/`. The bootstrap seeds are fixed in the code, so a rerun on the same inputs gives the same tables. Tests: `python -m unittest discover -s tests`.

## Licence

Code: MIT (see `LICENSE`). Third-party data and tool scores keep their own terms; files that may not be redistributed are not included, and the download list says how to obtain them.
