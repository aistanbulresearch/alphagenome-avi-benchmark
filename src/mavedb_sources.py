"""MaveDB part of the development positions (see data/development_sources.md).

For each saturation genome editing study listed in the source table, every score set of every listed
experiment is used, including all earlier versions reached through `supersededScoreSet`. From each
score set only two things are read: the `accession` column of its scores file (which variants were
scored; no score is read or kept) and the GRCh38 locations of its mapped variants. A scored variant
without a GRCh38 location is counted as unpositioned; the caller then applies the gene-region rule.

Usage: python src/mavedb_sources.py fetch     (records every file through src/download_public.py)
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from typing import Callable, Optional

from paths import ROOT
from paths import RAW as RAW_DIR
from clinvar_index import InputFormatError

API = "https://api.mavedb.org/api/v1"


class MissingSource(Exception):
    """A MaveDB record that was never obtained; the caller applies the gene-region rule."""


RAW = RAW_DIR
_LETTERS = "abcdefghijklmnopqrstuvwxyz"


def _exps(base: str, last: Optional[str] = None, meta: bool = True) -> tuple:
    letters = _LETTERS[:_LETTERS.index(last) + 1] if last else ""
    return tuple(([f"urn:mavedb:{base}-0"] if meta else []) + [f"urn:mavedb:{base}-{c}" for c in letters])


# Transcribed from data/development_sources.md (experiments=None: the gene region is used).
STUDIES = (
    {"study": "Findlay 2018", "gene": "BRCA1", "experiments": _exps("00000097", "z")},
    {"study": "Huang 2025", "gene": "BRCA2",
     "experiments": ("urn:mavedb:00001225-a", "urn:mavedb:00001279-a")},
    {"study": "Waters 2024", "gene": "BAP1", "experiments": _exps("00000662", "r")},
    {"study": "Olvera-Leon 2024", "gene": "RAD51C", "experiments": _exps("00000673", "i")},
    {"study": "Buckley 2024", "gene": "VHL", "experiments": ("urn:mavedb:00000675-a",)},
    {"study": "Radford 2023", "gene": "DDX3X", "experiments": _exps("00000658", "q") + ("urn:mavedb:00001228-a",)},
    {"study": "De Jonghe 2026", "gene": "RNU4-2", "experiments": None},
    {"study": "Lee 2025", "gene": "ATM", "experiments": None},
    {"study": "Dace 2025", "gene": "BRCA1", "experiments": None},
    {"study": "Woo 2025", "gene": "BARD1", "experiments": ("urn:mavedb:00001250-a",)},
    {"study": "Tejura 2026 PALB2", "gene": "PALB2", "experiments": ("urn:mavedb:00001259-a",)},
    {"study": "Tejura 2026 XRCC2", "gene": "XRCC2", "experiments": ("urn:mavedb:00001264-a",)},
    {"study": "Tejura 2026 SFPQ", "gene": "SFPQ", "experiments": ("urn:mavedb:00001265-a",)},
    {"study": "Tejura 2026 CTCF", "gene": "CTCF", "experiments": ("urn:mavedb:00001262-a",)},
)
_SUFFIX = {"experiment": "experiment.json", "score-set": "scoreset.json", "scores": "scores.csv",
           "mapped": "mapped.json"}


def local_name(urn: str, kind: str) -> str:
    """Path inside data/raw for one MaveDB record."""
    return f"mavedb/{urn.replace(':', '_')}.{_SUFFIX[kind]}"


def url(urn: str, kind: str) -> str:
    base = f"{API}/experiments/{urn}" if kind == "experiment" else f"{API}/score-sets/{urn}"
    return base + {"experiment": "", "score-set": "", "scores": "/scores", "mapped": "/mapped-variants"}[kind]


def _urn(value) -> Optional[str]:
    return value.get("urn") if isinstance(value, dict) else value


def score_set_urns(experiment: str, load: Callable[[str, str], dict]) -> list:
    """Every score set of an experiment: the listed ones and all versions they supersede."""
    found, todo = set(), list(load("experiment", experiment).get("scoreSetUrns") or [])
    while todo:
        urn = todo.pop()
        if urn in found:
            continue
        found.add(urn)
        older = _urn(load("score-set", urn).get("supersededScoreSet"))
        if older:
            todo.append(older)
    if not found:
        raise InputFormatError(f"{experiment}: no score set")
    return sorted(found)


def interval(location_holder: dict, chroms: dict) -> Optional[tuple]:
    """1-based inclusive GRCh38 span of a VRS Allele (interbase coordinates); None if not on GRCh38."""
    loc = (location_holder or {}).get("location") or {}
    chrom = chroms.get(((loc.get("sequenceReference") or {}).get("label")))
    start, end = loc.get("start"), loc.get("end")
    if chrom is None or not isinstance(start, int) or not isinstance(end, int):
        return None
    return (chrom, start, start + 1) if start == end else (chrom, start + 1, end)


def _alleles(post) -> list:
    """The Allele(s) of a post-mapped variant: itself, or the members of a haplotype or phased block."""
    if not isinstance(post, dict):
        return []
    if post.get("type") == "Allele":
        return [post]
    return post.get("members") or []


def score_set_positions(scores: Path, mapped: Path, chroms: dict) -> dict:
    """Intervals of every scored variant and the number of scored variants without a GRCh38 location."""
    with open(scores, encoding="utf-8", newline="") as handle:
        rows = csv.reader(line for line in handle if not line.startswith("#"))
        header = next(rows, [])
        if header.count("accession") != 1:
            raise InputFormatError(f"{scores.name}: no single 'accession' column")
        col = header.index("accession")
        scored, unusable = set(), 0
        for row in rows:
            accession = row[col].strip() if len(row) > col else ""
            if accession:
                scored.add(accession)
            else:
                unusable += 1                 # a scored row that cannot be identified, hence not located
    spans_of, complete, forms = {}, set(), {}
    for entry in json.loads(mapped.read_text(encoding="utf-8")):
        post = entry.get("postMapped")
        form = post.get("type", "untyped") if isinstance(post, dict) else "none"
        forms[form] = forms.get(form, 0) + 1
        alleles = _alleles(post)
        spans = [interval(a, chroms) for a in alleles]
        urn = entry.get("variantUrn")
        spans_of.setdefault(urn, set()).update(s for s in spans if s)
        if alleles and all(spans):
            complete.add(urn)                 # every member of the variant has a GRCh38 location
    intervals = {s for urn in scored for s in spans_of.get(urn, ())}
    unpositioned = unusable + sum(1 for urn in scored if urn not in complete)
    return {"scored": len(scored) + unusable, "unpositioned": unpositioned, "intervals": sorted(intervals),
            "forms": forms}


def _load_local(kind: str, urn: str) -> dict:
    path = RAW / local_name(urn, kind)
    if not path.exists():
        raise MissingSource(local_name(urn, kind))
    return json.loads(path.read_text(encoding="utf-8"))


def with_retry(download: Callable[[str, str], None], sleep: Callable[[float], None],
               attempts: int = 8, pause: float = 0.0) -> Callable[[str, str], None]:
    """Wrap a downloader: wait `pause` seconds before each request; on a transient failure (HTTP 429,
    502, 503 or 504, a time-out or a dropped connection) wait for the server's Retry-After (else 2, 4,
    8, ... seconds) and try again, at most `attempts` times. Other HTTP errors (404) are raised at once."""
    import urllib.error

    def wrapped(address: str, name: str) -> None:
        for attempt in range(attempts):
            if pause:
                sleep(pause)
            try:
                return download(address, name)
            except urllib.error.HTTPError as error:
                if error.code not in (429, 502, 503, 504) or attempt == attempts - 1:
                    raise
                after = error.headers.get("Retry-After") if error.headers is not None else None
                sleep(int(after) if after and str(after).isdigit() else 2 ** (attempt + 1))
            except (urllib.error.URLError, TimeoutError, ConnectionError):
                if attempt == attempts - 1:
                    raise
                sleep(2 ** (attempt + 1))
    return wrapped


def fetch(download: Optional[Callable[[str, str], None]] = None) -> list:
    """Record every needed MaveDB file (skipping files already recorded); return the records that
    MaveDB does not serve (HTTP 404), for which the builder applies the gene-region rule. Any other
    HTTP error stops the fetch."""
    import urllib.error
    if download is None:
        import time
        import download_public as dp

        def plain(address: str, name: str) -> None:
            dp.main(["--name", name, address])
        download = with_retry(plain, sleep=time.sleep, pause=1.0)

    unavailable: list = []

    def get(urn: str, kind: str) -> bool:
        name = local_name(urn, kind)
        if (RAW / name).exists():
            return True
        try:
            download(url(urn, kind), name)
        except urllib.error.HTTPError as error:
            if error.code != 404:
                raise
            unavailable.append(name)
            return False
        return True

    for study in STUDIES:
        for experiment in study["experiments"] or ():
            if not get(experiment, "experiment"):
                continue
            todo, seen = list(_load_local("experiment", experiment).get("scoreSetUrns") or []), set()
            while todo:
                urn = todo.pop()
                if urn in seen or not get(urn, "score-set"):
                    continue
                seen.add(urn)
                older = _urn(_load_local("score-set", urn).get("supersededScoreSet"))
                if older:
                    todo.append(older)
                get(urn, "scores")
                get(urn, "mapped")
    return unavailable


def study_score_sets(study: dict) -> list:
    """Score sets of a study from the recorded files (after fetch)."""
    return sorted({u for e in study["experiments"] or () for u in score_set_urns(e, _load_local)})


if __name__ == "__main__":
    if sys.argv[1:] == ["fetch"]:
        print(json.dumps({"unavailable": fetch()}, indent=1))
        sys.exit(0)
    print(__doc__)
    sys.exit(2)
