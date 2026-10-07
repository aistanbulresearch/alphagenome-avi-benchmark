"""Flags from ClinVar submission texts: does any submission of a variant mention a splicing predictor?

New intronic classifications often cite splicing predictions (for example SpliceAI under the ClinGen
splicing recommendations). Labels that rest on a predictor favour that predictor in a comparison. For each
variant of the variant list this script reads the free-text fields `Description` and
`ExplanationOfInterpretation` of every submission in submission_summary.txt.gz and records only whether
fixed terms occur. No number or other text is extracted or stored.

Output: data/derived/variants/submission_flags.tsv.gz with variation_id, submissions,
mentions_spliceai, mentions_any_splice_predictor (yes/no).

Usage:
    python src/submission_flags.py --submission-summary <submission_summary.txt.gz>
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

from paths import DERIVED

VARIANTS = DERIVED / "variants" / "clinvar_variants.tsv.gz"
OUT = DERIVED / "variants" / "submission_flags.tsv.gz"
SPLICEAI = re.compile(r"splice\s*-?\s*ai", re.IGNORECASE)
PREDICTORS = re.compile(
    r"splice\s*-?\s*ai|pangolin|maxent|max\s*ent\s*scan|spliceator|spip\b|splicesitefinder|nnsplice|"
    r"genesplicer|human splicing finder|\bhsf\b|dbscsnv|mmsplice|absplice|squirls|cryptic splice|"
    r"splic\w* predict\w*|predict\w* (?:to )?(?:\w+ ){0,3}splic\w*|in[\s-]silico\w*(?: \w+){0,4} splic\w*",
    re.IGNORECASE)


def flags(texts: list) -> tuple:
    joined = " ".join(texts)
    return bool(SPLICEAI.search(joined)), bool(PREDICTORS.search(joined))


def read_submissions(path: Path, wanted: set):
    """Yield (variation_id, description, explanation) for submissions of wanted variants."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as h:
        header = None
        for line in h:
            if line.startswith("#VariationID\t"):
                header = line[1:].rstrip("\n").split("\t")
                i_vid, i_desc, i_expl = (header.index(c) for c in ("VariationID", "Description",
                                                                    "ExplanationOfInterpretation"))
                continue
            if header is None or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) == len(header) and f[i_vid] in wanted:
                yield f[i_vid], f[i_desc], f[i_expl]
        if header is None:
            raise ValueError(f"{path.name}: no '#VariationID' header line")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--submission-summary", required=True, type=Path)
    a = ap.parse_args(argv)
    with gzip.open(VARIANTS, "rt", encoding="utf-8", newline="") as h:
        wanted = {r["variation_id"] for r in csv.DictReader(h, delimiter="\t")
                  if r["family"] in ("intronic", "utr") or r["region"] == "synonymous"}
    texts = defaultdict(list)
    for vid, desc, expl in read_submissions(a.submission_summary, wanted):
        texts[vid] += [desc, expl]
    with gzip.open(OUT, "wt", encoding="utf-8", newline="") as h:
        w = csv.writer(h, delimiter="\t", lineterminator="\n")
        w.writerow(("variation_id", "submissions", "mentions_spliceai", "mentions_any_splice_predictor"))
        for vid in sorted(wanted, key=int):
            s, p = flags(texts.get(vid, []))
            w.writerow((vid, len(texts.get(vid, [])) // 2, "yes" if s else "no", "yes" if p else "no"))
    print(f"variants {len(wanted)}; with submissions {len(texts)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
