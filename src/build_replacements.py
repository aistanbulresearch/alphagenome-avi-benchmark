"""ClinVar replacement history from the VCV XML release matching T1.

Streams the release once and reads, for each VariationArchive, only three direct children:
RecordStatus, ReplacedBy and ReplacedList/Replaced. Nested lists (for example RCV-level ReplacedList)
and every other element, including classifications, are ignored and discarded record by record.
Writes replaced_variation_id -> current_variation_id pairs (chains are followed by the counting
reader). A replaced identifier claimed by two different records, a replaced record without ReplacedBy,
an entry without a usable VCV identifier, or an accession that disagrees with its VariationID stops
the build. Never overwrites an existing output.

Usage: python src/build_replacements.py
"""
from __future__ import annotations

import gzip
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

from paths import ROOT
from paths import RAW as RAW_DIR
from clinvar_index import InputFormatError

RELEASE = RAW_DIR / "ClinVarVCVRelease_2026-10.xml.gz"
OUTPUT = ROOT / "data" / "derived" / "exposure" / "clinvar_replacements.tsv"
_VCV = re.compile(r"^VCV0*(\d+)$")


def _variation_id(entry: ET.Element) -> str:
    """The VariationID of a ReplacedBy or Replaced element, checked against its VCV accession."""
    attr = entry.get("VariationID")
    m = _VCV.match(entry.get("Accession", ""))
    from_accession = m.group(1) if m else None
    if attr is not None and not attr.isdigit():
        raise InputFormatError(f"non-numeric VariationID {attr!r}")
    if attr is not None and from_accession is not None and int(attr) != int(from_accession):
        raise InputFormatError(f"accession {entry.get('Accession')} disagrees with VariationID {attr}")
    vid = attr if attr is not None else from_accession
    if vid is None:
        raise InputFormatError(f"replacement entry without a VCV identifier ({entry.get('Accession')!r})")
    return str(int(vid))


def _pairs(archive: ET.Element):
    vid = archive.get("VariationID", "")
    if not vid.isdigit():
        raise InputFormatError("VariationArchive without a numeric VariationID")
    status = archive.findtext("RecordStatus", default="").strip()
    replaced_by = archive.find("ReplacedBy")
    if status == "replaced" and replaced_by is None:
        raise InputFormatError(f"record {vid} is replaced but names no ReplacedBy")
    if replaced_by is not None:
        yield vid, _variation_id(replaced_by)
    for entry in archive.findall("ReplacedList/Replaced"):       # direct child list only
        yield _variation_id(entry), vid


def build(release: Path = RELEASE, output: Path = OUTPUT) -> dict:
    """Write the pairs; return the number of records read and of distinct pairs."""
    if output.exists():
        raise FileExistsError(f"{output} exists; refusing to overwrite")
    mapping, records, root = {}, 0, None
    with gzip.open(release, "rb") as handle:
        for event, elem in ET.iterparse(handle, events=("start", "end")):
            if root is None:
                root = elem                                   # the release element
            if event == "end" and elem.tag == "VariationArchive":
                records += 1
                for old, new in _pairs(elem):
                    if old != new and mapping.setdefault(old, new) != new:
                        raise InputFormatError(f"VariationID {old} replaced by both {mapping[old]} and {new}")
                root.clear()                                  # discard the finished record entirely
    if not records:
        raise InputFormatError(f"{release.name}: no VariationArchive record")
    lines = ["replaced_variation_id\tcurrent_variation_id"]
    lines += [f"{old}\t{mapping[old]}" for old in sorted(mapping, key=int)]
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "x", encoding="utf-8", newline="\n") as out:
        out.write("\n".join(lines) + "\n")
    return {"records": records, "pairs": len(mapping)}


if __name__ == "__main__":
    print(build())
    sys.exit(0)
