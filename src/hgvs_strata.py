"""Deterministic stratum assignment from the HGVS c. expression in the ClinVar `Name` field.

Rules:

- Intronic family, by the signed offset from the exon boundary of the c. anchor (donor side '+',
  acceptor side '-'): canonical (+1, +2, -1, -2); donor_region (+3..+6); acceptor_region (-3..-20);
  near_intronic (+7..+50, -21..-50); deep_intronic (beyond +50 or -50). Offsets anchored in UTR
  positions (c.-N+K, c.*N+K) follow the same rules and are flagged utr_anchored.
- UTR family: utr5 (c.-N) and utr3 (c.*N) without an intronic offset.
- Coding family, from the protein expression: frameshift, extension, protein-level in-frame indel,
  then single-residue substitutions (synonymous, nonsense, missense); anything else coding_other;
  missing or uncertain protein consequence coding_unknown. Every coding-family record carries the
  flag annotation_pending (exonic splice-boundary status is decided separately).
- Ranges (deletions, duplications, inversions; insertions by their two flanking bases) use every
  affected base with precedence canonical > donor_region = acceptor_region > near_intronic >
  deep_intronic > coding = UTR. A range that leaves an intron or an exon crosses an exon-intron
  junction and is canonical. A tie between different strata at the highest precedence is
  annotation_ambiguous.
- Anything else (no single c. expression, n./g./m. expressions, multiple or uncertain expressions,
  unrecognised change) is unparsed, with a reason.

The parser reads only the `Name` string; it never uses positions from other fields.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Optional

PRECEDENCE = {"canonical": 5, "donor_region": 4, "acceptor_region": 4, "near_intronic": 3,
              "deep_intronic": 2, "coding": 1, "utr5": 1, "utr3": 1}
CODING_STRATA = ("frameshift", "extension", "coding_indel", "synonymous", "nonsense", "missense",
                 "coding_other", "coding_unknown")

_NAME = re.compile(r"^N[MR]_\d+(?:\.\d+)?(?:\([^)]*\))?:c\.(?P<body>\S+)(?:\s+\((?P<prot>p\.[^)]*\)?)\))?\s*$")
_POS = re.compile(r"^(?P<prefix>[-*]?)(?P<base>\d+)(?:(?P<sign>[+-])(?P<off>\d+))?")
_SUB = re.compile(r"^[ACGT]>[ACGT]$")
_RANGE_CHANGE = re.compile(r"^(?:del[ACGTN]*|dup[ACGTN]*|inv[ACGTN]*|ins[ACGTN]+|del[ACGTN]*ins[ACGTN]+)$")
_SINGLE_CHANGE = re.compile(r"^(?:del[ACGTN]?|dup[ACGTN]?|del[ACGTN]?ins[ACGTN]+)$")
_AA = "(?:Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val|Sec|Pyl)"
_SYN = re.compile(rf"^(?:(?:{_AA}|Ter)\d+)?=$")
_NONSENSE = re.compile(rf"^{_AA}\d+(?:Ter|\*)$")
_MISSENSE = re.compile(rf"^({_AA})\d+({_AA})$")


@dataclass(frozen=True)
class Position:
    """One c. position: prefix '' (coding), '-' (5' UTR) or '*' (3' UTR); optional intronic offset."""

    prefix: str
    base: int
    sign: Optional[str] = None
    offset: Optional[int] = None

    @property
    def intronic(self) -> bool:
        return self.offset is not None

    def stratum(self) -> str:
        if self.intronic:
            return intronic_stratum(self.sign, self.offset)
        return {"": "coding", "-": "utr5", "*": "utr3"}[self.prefix]

    def signed_offset(self) -> Optional[int]:
        if not self.intronic:
            return None
        return self.offset if self.sign == "+" else -self.offset


@dataclass(frozen=True)
class Assignment:
    """Family, stratum and the descriptive fields kept for reconstruction of other definitions."""

    family: str                          # intronic | utr | coding | unparsed | annotation_ambiguous
    stratum: str
    signed_offset: Optional[int] = None  # intronic only: closest-to-exon affected offset, signed
    utr_anchored: bool = False
    is_range: bool = False
    reason: str = ""
    flags: tuple = ()
    far_offset: Optional[int] = None      # intronic only: affected offset farthest from the exon, signed;
                                          # None when a range spans the middle of an intron
    ends: Optional[tuple] = None          # intronic ranges: signed offsets of both ends (None for an exonic end)


def intronic_stratum(sign: str, offset: int) -> str:
    """Stratum of one intronic base given its side ('+' donor, '-' acceptor) and offset (>= 1)."""
    if offset < 1:
        raise ValueError("intronic offset must be at least 1")
    if offset <= 2:
        return "canonical"
    if sign == "+":
        if offset <= 6:
            return "donor_region"
        return "near_intronic" if offset <= 50 else "deep_intronic"
    if offset <= 20:
        return "acceptor_region"
    return "near_intronic" if offset <= 50 else "deep_intronic"


def parse_position(text: str) -> tuple[Optional[Position], str]:
    """Parse one c. position at the start of `text`; return it and the remaining text."""
    m = _POS.match(text)
    if not m or int(m.group("base")) == 0:
        return None, text
    off = m.group("off")
    if off is not None and int(off) == 0:
        return None, text
    pos = Position(m.group("prefix"), int(m.group("base")), m.group("sign"), int(off) if off else None)
    return pos, text[m.end():]


def _next_exon_start(pos: Position) -> tuple[str, int]:
    """Anchor (prefix, base) of the exon that follows the intron after a donor-side anchor."""
    if pos.prefix == "-":
        return ("", 1) if pos.base == 1 else ("-", pos.base - 1)
    return pos.prefix, pos.base + 1


def _affected(start: Position, end: Position) -> Optional[list[Position]]:
    """Affected positions that can decide the stratum of a range; None if they cannot be resolved."""
    if not start.intronic and not end.intronic:
        return [start, end]
    if start.intronic and end.intronic:
        if (start.prefix, start.base, start.sign) == (end.prefix, end.base, end.sign):
            lo, hi = sorted((start.offset, end.offset))
            return [Position(start.prefix, start.base, start.sign, k) for k in range(lo, hi + 1)]
        if start.sign == "+" and end.sign == "-" and _next_exon_start(start) == (end.prefix, end.base):
            # Within one intron: the bases closest to either exon decide; deeper bases rank lower.
            return [start, end]
    # The range leaves an intron or an exon, so it contains offset 1 of the first junction it meets:
    # an acceptor if it starts inside an intron, a donor if it starts inside an exon.
    sign = "-" if start.intronic else "+"
    return [start, end, Position(start.prefix, start.base, sign, 1)]


def _decide(bases: list[Position]) -> Optional[Position]:
    """Highest-precedence base; None if bases of different strata tie at the top."""
    top = max(PRECEDENCE[p.stratum()] for p in bases)
    leaders = [p for p in bases if PRECEDENCE[p.stratum()] == top]
    if len({p.stratum() for p in leaders}) > 1:
        return None
    leaders.sort(key=lambda p: (p.offset or 0, p.sign != "+"))
    return leaders[0]


def protein_class(prot: Optional[str]) -> str:
    """Coding-family stratum from the protein expression, given with or without the leading 'p.'."""
    if not prot:
        return "coding_unknown"
    p = prot.strip()
    if p.startswith("p."):
        p = p[2:]
    if p.startswith("(") and p.endswith(")"):
        p = p[1:-1]
    if not p:
        return "coding_unknown"
    if "fs" in p:
        return "frameshift"
    if "ext" in p:
        return "extension"
    if "?" in p:
        return "coding_unknown"
    if any(tag in p for tag in ("del", "ins", "dup")):
        return "coding_indel"
    if _SYN.match(p):
        return "synonymous"
    if _NONSENSE.match(p):
        return "nonsense"
    m = _MISSENSE.match(p)
    if m and m.group(1) != m.group(2):
        return "missense"
    return "coding_other"


def _unparsed(reason: str) -> Assignment:
    return Assignment("unparsed", "unparsed", reason=reason)


def assign(name: Optional[str]) -> Assignment:
    """Assign one ClinVar `Name` value to a family and stratum."""
    if not name or not name.strip():
        return _unparsed("missing name")
    m = _NAME.match(name.strip())
    if not m:
        return _unparsed("no single c. expression on an NM_ or NR_ transcript")
    body = m.group("body")
    if body.startswith(("[", "(")) or "?" in body or ";" in body:
        return _unparsed("multiple or uncertain expression")
    start, rest = parse_position(body)
    if start is None:
        return _unparsed("no position")
    end, is_range = start, rest.startswith("_")
    if is_range:
        end, rest = parse_position(rest[1:])
        if end is None:
            return _unparsed("unresolved range end")
        if not _RANGE_CHANGE.match(rest):
            return _unparsed("unrecognised range change")
        if _last_intron_unknown(start, end):
            return Assignment("annotation_ambiguous", "annotation_ambiguous", is_range=True,
                              reason="range may lie within the intron before the 3' UTR")
    elif not (_SUB.match(rest) or _SINGLE_CHANGE.match(rest)):
        return _unparsed("unrecognised change")
    bases = _affected(start, end) if is_range else [start]
    chosen = _decide(bases)
    if chosen is None:
        return Assignment("annotation_ambiguous", "annotation_ambiguous", is_range=is_range,
                          reason="different strata tie at the highest precedence")
    utr_anchored = any(p.intronic and p.prefix in ("-", "*") for p in bases)
    stratum = chosen.stratum()
    if chosen.intronic:
        return Assignment("intronic", stratum, chosen.signed_offset(), utr_anchored, is_range,
                          far_offset=_far_offset(start, end, bases, is_range),
                          ends=(start.signed_offset(), end.signed_offset()) if is_range else None)
    if stratum in ("utr5", "utr3"):
        return Assignment("utr", stratum, is_range=is_range)
    return Assignment("coding", protein_class(m.group("prot")), is_range=is_range,
                      flags=("annotation_pending",))


def _far_offset(start: Position, end: Position, bases: list[Position], is_range: bool) -> Optional[int]:
    """Signed offset of the affected intronic base farthest from the exon (see Assignment.far_offset)."""
    if is_range and start.intronic and end.intronic and start.sign == "+" and end.sign == "-" \
            and _next_exon_start(start) == (end.prefix, end.base):
        return None
    intronic = [p for p in bases if p.intronic]
    return max(intronic, key=lambda p: p.offset).signed_offset()


def _last_intron_unknown(start: Position, end: Position) -> bool:
    """True for c.N+a_*1-b: whether N is the last coding base cannot be told from the name."""
    return (start.prefix == "" and start.sign == "+" and end.prefix == "*" and end.base == 1
            and end.sign == "-")


__all__ = ["Assignment", "CODING_STRATA", "PRECEDENCE", "Position", "assign", "intronic_stratum",
           "parse_position", "protein_class"]
