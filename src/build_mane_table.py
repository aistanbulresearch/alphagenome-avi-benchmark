"""MANE Select transcript table for the binding strata.

Reads the MANE release 1.5 summary (a recorded preparation file) and writes, for every MANE Select
transcript, its RefSeq accession.version, its NCBI GeneID and whether the gene also has a MANE Plus
Clinical transcript. Fails closed on any unexpected row. Never overwrites an existing output.

Usage: python src/build_mane_table.py
"""
from __future__ import annotations

import gzip
import re
import sys
from collections import Counter
from pathlib import Path

from paths import ROOT
from paths import RAW as RAW_DIR
from clinvar_index import InputFormatError

SUMMARY = RAW_DIR / "MANE.GRCh38.v1.5.summary.txt.gz"
OUTPUT = ROOT / "data" / "derived" / "exposure" / "mane_select.tsv"
REQUIRED = ("#NCBI_GeneID", "RefSeq_nuc", "MANE_status")
STATUSES = ("MANE Select", "MANE Plus Clinical")
_TRANSCRIPT = re.compile(r"^N[MR]_\d+\.\d+$")
_GENE = re.compile(r"^GeneID:(\d+)$")


def _rows(summary: Path):
    with gzip.open(summary, "rt", encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        missing = [f for f in REQUIRED if header.count(f) != 1]
        if missing:
            raise InputFormatError(f"{summary.name}: header lacks or repeats {', '.join(missing)}")
        cols = {f: header.index(f) for f in REQUIRED}
        for line_no, line in enumerate(handle, 2):
            values = line.rstrip("\r\n").split("\t")
            if len(values) != len(header):
                raise InputFormatError(f"{summary.name}: malformed row at line {line_no}")
            gene, transcript, status = (values[cols[f]] for f in REQUIRED)
            if status not in STATUSES or not _TRANSCRIPT.match(transcript) or not _GENE.match(gene):
                raise InputFormatError(f"{summary.name}: unexpected value at line {line_no}")
            yield _GENE.match(gene).group(1), transcript, status


def build(summary: Path = SUMMARY, output: Path = OUTPUT) -> dict:
    """Write the table; return the row counts by status."""
    if output.exists():
        raise FileExistsError(f"{output} exists; refusing to overwrite")
    select, plus, counts = {}, set(), Counter()
    for gene, transcript, status in _rows(summary):
        counts[status] += 1
        if status == "MANE Plus Clinical":
            plus.add(gene)
        elif gene in select:
            raise InputFormatError(f"{summary.name}: gene {gene} has two MANE Select transcripts")
        else:
            select[gene] = transcript
    if not select:
        raise InputFormatError(f"{summary.name}: no MANE Select row")
    if plus - set(select):
        raise InputFormatError(f"{summary.name}: a MANE Plus Clinical gene has no MANE Select transcript")
    lines = ["refseq_transcript\tgene_id\tplus_clinical"]
    lines += [f"{tx}\t{gene}\t{'yes' if gene in plus else 'no'}"
              for gene, tx in sorted(select.items(), key=lambda kv: kv[1])]
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "x", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(lines) + "\n")
    return {**dict(counts), "genes with Plus Clinical": len(plus)}


if __name__ == "__main__":
    print(build())
    sys.exit(0)
