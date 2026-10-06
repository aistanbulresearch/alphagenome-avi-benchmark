"""Small derived exposure files.

- development_submitters.txt: the listed organization (checked against ClinVar's organization
  list, a preparation file) and the listed spellings of its name.
- hg19ToHg38.over.chain.gz: the UCSC chain, copied byte for byte from its recorded download.
- development_sources.tsv: each development source with its single permitted status and its evidence
  file; written last, after the positions, and checked with the coverage rule.
Outputs are never overwritten.

Usage: python src/build_exposure_files.py submitters | chain
       (the sources file is written by build_development_positions.py)
"""
from __future__ import annotations

import hashlib
import shutil
import sys
from pathlib import Path

from paths import ROOT, sha256_file
from paths import RAW as RAW_DIR
from development_sources import (REQUIRED_DEV_SOURCES, SOURCE_STATUSES, SUBMITTER_ORGANIZATIONS,
                               SUBMITTER_SPELLINGS, coverage_problems)
from clinvar_index import InputFormatError

RAW = RAW_DIR
EXPOSURE = ROOT / "data" / "derived" / "exposure"
ORGANIZATIONS = RAW / "organization_summary.txt"
CHAIN = RAW / "hg19ToHg38.over.chain.gz"


def _organization_names(path: Path) -> dict:
    """Organization ID -> list of names in ClinVar's organization list."""
    out: dict = {}
    with open(path, encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        if header[:2] != ["#organization", "organization ID"]:
            raise InputFormatError(f"{path.name}: unexpected header")
        for line in handle:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) >= 2:
                out.setdefault(parts[1].strip(), []).append(parts[0].strip())
    return out


def write_submitters(organizations: Path = ORGANIZATIONS,
                     output: Path = EXPOSURE / "development_submitters.txt") -> None:
    names = _organization_names(organizations)
    for org, name in SUBMITTER_ORGANIZATIONS.items():
        if names.get(org) != [name]:
            raise InputFormatError(f"organization {org} is not listed once under its listed name")
    lines = ["organization_id\tsubmitter_name"]
    lines += [f"{org}\t{spelling}" for org in sorted(SUBMITTER_SPELLINGS) for spelling in sorted(SUBMITTER_SPELLINGS[org])]
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "x", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(lines) + "\n")


def copy_chain(source: Path = CHAIN, output: Path = EXPOSURE / "hg19ToHg38.over.chain.gz") -> str:
    """Copy the chain unchanged; return its SHA-256 (equal on both sides)."""
    if output.exists():
        raise FileExistsError(f"{output} exists; refusing to overwrite")
    output.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, output)
    digest = sha256_file(source)
    if hashlib.sha256(output.read_bytes()).hexdigest() != digest:
        output.unlink()
        raise IOError("chain copy differs from its source")
    return digest


def write_sources(evidence: dict, output: Path = EXPOSURE / "development_sources.tsv") -> None:
    """Write source, status, evidence for every required source; refuse unless the
    coverage rule passes on the result."""
    if set(evidence) != set(REQUIRED_DEV_SOURCES):
        raise InputFormatError("evidence must be given for exactly the required development sources")
    rows = {s: (SOURCE_STATUSES[s][0], evidence[s]) for s in REQUIRED_DEV_SOURCES}
    problems = coverage_problems(rows, output.parent)
    if problems:
        raise InputFormatError("; ".join(problems))
    lines = ["source\tstatus\tevidence"] + [f"{s}\t{rows[s][0]}\t{rows[s][1]}" for s in REQUIRED_DEV_SOURCES]
    with open(output, "x", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    task = sys.argv[1] if len(sys.argv) > 1 else ""
    if task == "submitters":
        write_submitters()
    elif task == "chain":
        print(copy_chain())
    else:
        print(__doc__)
        sys.exit(2)
    sys.exit(0)
