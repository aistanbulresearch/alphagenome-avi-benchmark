"""The developers' benchmark and training sources whose variants are excluded from the new-label set.

    complete         the source's variant list was obtained; its positions are excluded
    submitter rule   GREGoR has no published variant list: ClinVar variants with a submission from the
                     organisation that supplied the developers' GREGoR cases are excluded
    not excluded     evaluated by the developers only on a different task (trait association)
"""
from __future__ import annotations

from pathlib import Path

from clinvar_index import InputFormatError

SOURCE_STATUSES = {
    "saturation genome editing": ("complete",),
    "TraitGym Mendelian": ("complete",),
    "TraitGym Complex": ("complete",),
    "ncVarDB": ("complete",),
    "GREGoR": ("submitter rule",),
    "complex-trait benchmark": ("not excluded",),
}
REQUIRED_DEV_SOURCES = tuple(SOURCE_STATUSES)
# The organisation, and the spellings of its name found in ClinVar: the name in ClinVar's organisation list
# and the same name with ', ' in place of '; '.
SUBMITTER_ORGANIZATIONS = {"506627": "Broad Center for Mendelian Genomics; Broad Institute of MIT and Harvard"}
SUBMITTER_SPELLINGS = {org: {name, name.replace("; ", ", ")} for org, name in SUBMITTER_ORGANIZATIONS.items()}


def read_sources(path: Path) -> dict:
    """Source coverage table: TSV with header source, status, evidence."""
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0].split("\t")[:3] != ["source", "status", "evidence"]:
        raise InputFormatError(f"{path.name}: header must be source, status, evidence")
    out = {}
    for line in lines[1:]:
        if line.strip():
            parts = line.split("\t")
            out[parts[0]] = (parts[1] if len(parts) > 1 else "", parts[2] if len(parts) > 2 else "")
    return out


def coverage_problems(sources: dict, base: Path) -> list[str]:
    """Every source has one of its statuses and an existing evidence file inside `base`."""
    problems = []
    for name in REQUIRED_DEV_SOURCES:
        status, evidence = sources.get(name, ("", ""))
        target = (base / evidence).resolve() if evidence else None
        if status not in SOURCE_STATUSES[name]:
            problems.append(f"source without a permitted status ({', '.join(SOURCE_STATUSES[name])}): {name}")
        elif not evidence:
            problems.append(f"evidence file not named for source: {name}")
        elif Path(evidence).is_absolute() or base.resolve() not in target.parents:
            problems.append(f"evidence file outside the exposure folder: {name}")
        elif not target.is_file():
            problems.append(f"evidence file missing for source: {name}")
    return problems


__all__ = ["SOURCE_STATUSES", "REQUIRED_DEV_SOURCES", "SUBMITTER_ORGANIZATIONS", "SUBMITTER_SPELLINGS",
           "read_sources", "coverage_problems"]
