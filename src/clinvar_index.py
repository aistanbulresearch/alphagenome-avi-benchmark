"""T1 candidate assembly and the exposure index.

Streaming readers for gzipped variant_summary, VCF and submission_summary files; the T1 file is reduced
to one candidate per VariationID under the duplicate rule; every row of every exposure source is then
matched to the candidates by GRCh38 allele, GRCh37 allele, AlleleID or VariationID (each VariationID
mapped through ClinVar's replacement history) before any other filter. Exposure files fail closed: a malformed exposure row stops the run,
because a lost exposure row could make a seen variant look new. Standard library only.
"""
from __future__ import annotations

import bisect
import gzip
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Iterator, Optional

import clinvar_records as cr

MAX_UNKNOWN_SPAN = 50          # bases; see ExposureIndex.block_unmatched
BASE_EXPOSURE = ("E-dev", "T0", "E-hist")
VCF_COLUMNS = ("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
VCF_SINGLE_INFO = ("ALLELEID", "CLNSIG", "CLNSIGINCL", "CLNREVSTAT")


class InputFormatError(RuntimeError):
    """A required header field is missing or duplicated, or an exposure row is malformed."""


def _check_header(header: list[str], required: Iterable[str], path: Path) -> None:
    missing = [f for f in required if f not in header]
    duplicated = sorted({f for f in header if header.count(f) > 1})
    if missing or duplicated:
        raise InputFormatError(f"{path.name}: missing fields {missing}, duplicated fields {duplicated}")


def read_variant_summary(path: Path, source: str, flow: Counter, fail_closed: bool = True) -> Iterator[cr.Record]:
    """One Record per data row. A malformed row stops the run if `fail_closed`, else it is counted."""
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        _check_header(header, cr.REQUIRED_VS_FIELDS, path)
        for line_no, line in enumerate(handle, 2):
            values = line.rstrip("\r\n").split("\t")
            flow[f"{source}: rows read"] += 1
            if len(values) != len(header):
                if fail_closed:
                    raise InputFormatError(f"{source}: malformed row at line {line_no} (field count)")
                flow[f"{source}: malformed row (field count), excluded"] += 1
                continue
            yield cr.from_variant_summary(dict(zip(header, values)), source)


def check_info(info: str, source: str, line_no: int) -> None:
    """ALLELEID present once with a value; classification fields at most once, with valid syntax."""
    items = [item for item in info.split(";") if item]
    keys = [item.split("=", 1)[0] for item in items]
    for k in VCF_SINGLE_INFO:
        if keys.count(k) > 1:
            raise InputFormatError(f"{source}: VCF line {line_no} repeats INFO field {k}")
    values = {item.split("=", 1)[0]: (item.split("=", 1)[1] if "=" in item else None) for item in items}
    if not values.get("ALLELEID"):
        raise InputFormatError(f"{source}: VCF line {line_no} has no ALLELEID value")
    for k in ("CLNSIG", "CLNSIGINCL", "CLNREVSTAT"):
        if k in values and values[k] is None:
            raise InputFormatError(f"{source}: VCF line {line_no} has INFO field {k} without '='")
        if k in values and not values[k] and k != "CLNSIG":      # 'CLNSIG=' means no classification
            raise InputFormatError(f"{source}: VCF line {line_no} has INFO field {k} without a value")
    if values.get("CLNSIGINCL"):
        try:
            cr.incl_pairs(values["CLNSIGINCL"])
        except ValueError:
            raise InputFormatError(f"{source}: VCF line {line_no} has a malformed CLNSIGINCL pair") from None


def read_vcf(path: Path, source: str, flow: Counter) -> Iterator[cr.Record]:
    """Records from a GRCh38 ClinVar VCF; fails closed on a missing header or a malformed data line."""
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        width = None
        for line_no, line in enumerate(handle, 1):
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                columns = line.rstrip("\r\n").split("\t")
                if columns[:8] != list(VCF_COLUMNS):
                    raise InputFormatError(f"{path.name}: VCF header columns {columns[:8]}")
                width = len(columns)
                continue
            fields = line.rstrip("\r\n").split("\t")
            if width is None or len(fields) != width:
                raise InputFormatError(f"{source}: malformed VCF line {line_no}")
            check_info(fields[7], source, line_no)
            flow[f"{source}: rows read"] += 1
            yield from cr.from_vcf_line(line, source)
        if width is None:
            raise InputFormatError(f"{source}: no #CHROM header line")


def read_submitters(path: Path, flow: Counter) -> Iterator[tuple[str, str]]:
    """(VariationID, Submitter) pairs from a submission_summary file; fails closed."""
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        header = None
        for line_no, line in enumerate(handle, 1):
            if header is None:
                if line.startswith("#VariationID\t"):      # not the '#VariationID:' documentation line
                    header = line.rstrip("\r\n").split("\t")
                    _check_header(header, ("#VariationID", "Submitter"), path)
                    i_vid, i_sub = header.index("#VariationID"), header.index("Submitter")
                elif not line.startswith("#"):              # only documentation lines precede the header
                    raise InputFormatError(f"submission_summary: data before the header at line {line_no}")
                continue
            values = line.rstrip("\r\n").split("\t")
            if len(values) != len(header):
                raise InputFormatError(f"submission_summary: malformed row at line {line_no}")
            flow["submission_summary: rows read"] += 1
            yield values[i_vid].strip(), values[i_sub].strip()
        if header is None:
            raise InputFormatError("submission_summary: header line not found")


def read_replacements(path: Path) -> dict:
    """ClinVar replacement history: TSV with header replaced_variation_id, current_variation_id."""
    out = {}
    with open(path, encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        if header[:2] != ["replaced_variation_id", "current_variation_id"]:
            raise InputFormatError(f"{path.name}: header must be replaced_variation_id, current_variation_id")
        for line_no, line in enumerate(handle, 2):
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 2 or not parts[0] or not parts[1]:
                raise InputFormatError(f"{path.name}: malformed line {line_no}")
            out[parts[0]] = parts[1]
    # Follow chains (a -> b -> c) to the current identifier; a cycle is a format error.
    resolved = {}
    for start in out:
        seen, v = {start}, out[start]
        while v in out:
            if v in seen:
                raise InputFormatError(f"{path.name}: replacement cycle at {start}")
            seen.add(v)
            v = out[v]
        resolved[start] = v
    return resolved


@dataclass
class Candidate:
    """One T1 VariationID that passed the duplicate rule, with its exposure state."""

    record: cr.Record
    group: str
    keys37: frozenset = frozenset()
    spans37: frozenset = frozenset()
    unresolved: Optional[str] = None
    development: Optional[str] = None     # reason for exclusion as AVI development material
    present: set = field(default_factory=set)
    seen: set = field(default_factory=set)
    seen_strict: set = field(default_factory=set)     # sources that saw a strict term (not only unknown)
    max_exposure_stars: int = -1


@dataclass
class T1Universe:
    """Identifiers and alleles of every T1 row (missing values excluded; VariationIDs mapped to their
    current replacement), for the unknown-record diagnosis."""

    canon: Callable[[Optional[str]], Optional[str]] = lambda v: v
    variation_ids: set = field(default_factory=set)
    allele_ids: set = field(default_factory=set)
    keys38: set = field(default_factory=set)
    keys37: set = field(default_factory=set)

    def add(self, rec: cr.Record) -> None:
        for target, value in ((self.variation_ids, self.canon(rec.variation_id)), (self.allele_ids, rec.allele_id),
                              (self.keys38, rec.key), (self.keys37, rec.key37)):
            if value is not None:
                target.add(value)

    def knows(self, rec: cr.Record) -> bool:
        return any(v is not None and v in s for v, s in (
            (self.canon(rec.variation_id), self.variation_ids), (rec.allele_id, self.allele_ids),
            (rec.key, self.keys38), (rec.key37, self.keys37)))


def assemble_t1(records: Iterator[cr.Record], flow: Counter,
                universe: Optional[T1Universe] = None) -> dict[str, Candidate]:
    """Candidates strict at T1, one per VariationID; exclusions counted in `flow`."""
    strict_rows: dict[str, list[cr.Record]] = defaultdict(list)
    keys37: dict[str, set] = defaultdict(set)
    spans37: dict[str, set] = defaultdict(set)
    nonstrict_vids: set[str] = set()
    strict_any_assembly: set[str] = set()
    for rec in records:
        if universe is not None:
            universe.add(rec)
        group = cr.strict_group(rec.classification)
        if rec.variation_id is None:
            flow["T1: row without VariationID"] += 1
            continue
        if rec.key37 is not None:
            keys37[rec.variation_id].add(rec.key37)      # every GRCh37 row of the VariationID
        if rec.assembly == "GRCh37" and rec.span is not None:
            spans37[rec.variation_id].add(rec.span)
        if group is not None:
            strict_any_assembly.add(rec.variation_id)
        if rec.assembly != "GRCh38":
            continue
        if group is None:
            nonstrict_vids.add(rec.variation_id)
        else:
            strict_rows[rec.variation_id].append(rec)
    flow["T1 strict VariationIDs: no GRCh38 row"] = len(strict_any_assembly - set(strict_rows))
    candidates: dict[str, Candidate] = {}
    for vid, rows in strict_rows.items():
        if len({r.identity() for r in rows}) > 1 or vid in nonstrict_vids:
            flow["T1 strict VariationIDs: ambiguous_duplicate"] += 1
            continue
        rec = rows[0]
        if rec.key is None:
            flow["T1 strict VariationIDs: no GRCh38 VCF-style allele"] += 1
            continue
        if rec.allele_id is None:
            flow["T1 strict VariationIDs: no AlleleID"] += 1
            continue
        candidates[vid] = Candidate(rec, cr.strict_group(rec.classification), frozenset(keys37.get(vid, ())),
                                    frozenset(spans37.get(vid, ())))
    flow["T1 strict VariationIDs: candidates"] = len(candidates)
    _mark_internal_collisions(candidates)
    return candidates


def _mark_internal_collisions(candidates: dict[str, Candidate]) -> None:
    by_key, by_allele = defaultdict(set), defaultdict(set)
    for c in candidates.values():
        by_key[c.record.key].add(c.record.allele_id)
        by_allele[c.record.allele_id].add(c.record.key)
    for c in candidates.values():
        if len(by_key[c.record.key]) > 1 or len(by_allele[c.record.allele_id]) > 1:
            c.unresolved = "identity collision at T1"


def key_span(key: str) -> tuple[str, int, int]:
    chrom, pos, ref, _ = key.split(":")
    return chrom, int(pos), int(pos) + len(ref) - 1


class ExposureIndex:
    """Matches exposure records to candidates by any key, with every VariationID mapped
    to its current replacement on both sides; records collisions and the strict exposure records of
    which no row is known at T1 by any identifier or allele."""

    def __init__(self, candidates: dict[str, Candidate], universe: Optional[T1Universe] = None,
                 aliases: Optional[dict] = None):
        self.candidates = candidates
        self.universe = universe
        self.aliases = aliases or {}
        self.by_key, self.by_key37, self.by_allele, self.by_vid = (defaultdict(set) for _ in range(4))
        self.diagnosis: dict[str, dict] = defaultdict(dict)
        self._rows = 0
        for vid, c in candidates.items():
            self.by_key[c.record.key].add(vid)
            self.by_allele[c.record.allele_id].add(vid)
            self.by_vid[self.canon(vid)].add(vid)
            for k in c.keys37:
                self.by_key37[k].add(vid)

    def canon(self, vid: Optional[str]) -> Optional[str]:
        return None if vid is None else self.aliases.get(vid, vid)

    def add(self, rec: cr.Record, flow: Counter) -> None:
        matched = set()
        if rec.key is not None:
            matched |= self.by_key.get(rec.key, set())
        if rec.key37 is not None:
            matched |= self.by_key37.get(rec.key37, set())
        if rec.allele_id is not None:
            matched |= self.by_allele.get(rec.allele_id, set())
        if rec.variation_id is not None:
            matched |= self.by_vid.get(self.canon(rec.variation_id), set())
        for term in rec.unknown_terms:
            flow[f"{rec.source}: unknown classification term (counted as seen): {term}"] += 1
        if rec.seen and self.universe is not None:
            self._diagnose(rec)
        if not matched:
            return
        flow[f"{rec.source}: rows matching a T1 candidate"] += 1
        for vid in matched:
            c = self.candidates[vid]
            c.present.add(rec.source)
            if rec.seen:
                c.seen.add(rec.source)
                if rec.strict_token:
                    c.seen_strict.add(rec.source)
                if isinstance(rec.stars, int):
                    c.max_exposure_stars = max(c.max_exposure_stars, rec.stars)
            if rec.assembly == "GRCh38" and rec.key is not None and rec.allele_id is not None:
                if (rec.allele_id == c.record.allele_id) != (rec.key == c.record.key) and c.unresolved is None:
                    c.unresolved = f"identity collision with {rec.source}"

    def _diagnose(self, rec: cr.Record) -> None:
        """Aggregate, per source and identity, whether any row is known at T1 and where it lies. Records
        are grouped only by a real identifier; a record without one is its own group."""
        self._rows += 1
        ident = (self.canon(rec.variation_id) and f"vid:{self.canon(rec.variation_id)}") or \
                (rec.allele_id and f"aid:{rec.allele_id}") or (rec.key and f"key:{rec.key}") or \
                (rec.key37 and f"key37:{rec.key37}") or f"row:{rec.source}:{self._rows}"
        entry = self.diagnosis[rec.source].setdefault(ident, {"known": False, "pos38": set(), "pos37": set()})
        if self.universe.knows(rec):
            entry["known"] = True
        if rec.span is not None:
            entry["pos38" if rec.assembly == "GRCh38" else "pos37"].add(rec.span)
        if rec.key is not None:
            entry["pos38"].add(key_span(rec.key))
        if rec.key37 is not None:
            entry["pos37"].add(key_span(rec.key37))

    def unmatched(self) -> dict:
        """Identity -> positions for strict exposure records unknown at T1 (all sources)."""
        out = {}
        for source, entries in self.diagnosis.items():
            for ident, e in entries.items():
                if not e["known"]:
                    prev = out.setdefault(ident, {"pos38": set(), "pos37": set(), "sources": set()})
                    prev["pos38"] |= e["pos38"]
                    prev["pos37"] |= e["pos37"]
                    prev["sources"].add(source)
        return out

    def block_unmatched(self, lift: Optional[Callable] = None) -> dict:
        """Exclude from the new-label set every candidate at the position of an unknown strict exposure
        record of at most MAX_UNKNOWN_SPAN bases (a small variant that may be the same variant under an
        older identifier; a long structural record says nothing about a small variant inside it):
        GRCh38 positions directly; GRCh37 positions against the candidate's GRCh37 rows and, through
        `lift` (GRCh37 span -> list of GRCh38 spans), against its GRCh38 span. Returns counts."""
        unknown = self.unmatched()
        small = lambda p: p[2] - p[1] + 1 <= MAX_UNKNOWN_SPAN  # noqa: E731
        spans38 = [p for u in unknown.values() for p in u["pos38"] if small(p)]
        spans37 = [p for u in unknown.values() for p in u["pos37"] if small(p)]
        long_records = sum(1 for u in unknown.values()
                           if any(not small(p) for p in u["pos38"] | u["pos37"]))
        unlifted = partly = 0
        if lift is not None:
            for p in spans37:
                mapped = lift(p)
                spans38.extend(mapped)
                covered = sum(e - s + 1 for _, s, e in mapped)
                if not mapped:
                    unlifted += 1
                elif covered < p[2] - p[1] + 1:
                    partly += 1
        idx38, idx37 = PositionIndex(spans38), PositionIndex(spans37)
        reason = "unresolved historical identity (same position)"
        for c in self.candidates.values():
            hit = idx38.overlaps(*key_span(c.record.key))
            hit = hit or any(idx37.overlaps(*k) for k in [key_span(k) for k in c.keys37] + list(c.spans37))
            if hit and c.unresolved is None:
                c.unresolved = reason
        blocked = sum(1 for c in self.candidates.values() if c.unresolved == reason)
        placed = sum(1 for u in unknown.values() if u["pos38"] or u["pos37"])
        return {"unknown_records": len(unknown), "placed": placed, "unplaced": len(unknown) - placed,
                "records_with_a_span_over_limit": long_records,
                "grch37_positions_not_converted": unlifted if lift is not None else None,
                "grch37_positions_partly_converted": partly if lift is not None else None,
                "candidates_blocked": blocked}


class PositionIndex:
    """Merged intervals (1-based, inclusive) per chromosome."""

    def __init__(self, intervals: Iterable[tuple[str, int, int]]):
        per = defaultdict(list)
        self.starts, self.ends = {}, {}
        for chrom, start, end in intervals:
            per[cr.canonical_chrom(chrom)].append((min(start, end), max(start, end)))
        for chrom, ivs in per.items():
            ivs.sort()
            merged = [list(ivs[0])]
            for s, e in ivs[1:]:
                if s <= merged[-1][1] + 1:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])
            self.starts[chrom] = [s for s, _ in merged]
            self.ends[chrom] = [e for _, e in merged]

    def overlaps(self, chrom: str, start: int, end: int) -> bool:
        c = cr.canonical_chrom(chrom)
        if c not in self.starts:
            return False
        i = bisect.bisect_right(self.starts[c], end) - 1
        return i >= 0 and self.ends[c][i] >= start


def read_positions(path: Path) -> dict[str, list[tuple[str, int, int]]]:
    """Development-benchmark positions: TSV with header chrom, start, end, source (GRCh38, 1-based)."""
    out = defaultdict(list)
    with open(path, encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        if header[:4] != ["chrom", "start", "end", "source"]:
            raise InputFormatError(f"{path.name}: header must be chrom, start, end, source")
        for line_no, line in enumerate(handle, 2):
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 4:
                raise InputFormatError(f"{path.name}: malformed line {line_no}")
            out[parts[3]].append((parts[0], int(parts[1]), int(parts[2])))
    return out


def candidate_span(c: Candidate) -> tuple[str, int, int]:
    """GRCh38 bases spanned by the candidate's reference allele (an insertion spans its anchor base)."""
    return key_span(c.record.key)


def mark_development(candidates: dict[str, Candidate], positions: dict[str, list],
                     excluded_vids: dict[str, set], canon: Callable = lambda v: v) -> Counter:
    """Mark candidates that are AVI development material: overlap with a benchmark position (any
    allele), or a VariationID listed by a development-submitter source (both sides mapped to the current
    replacement by `canon`). Returns counts per source."""
    counts = Counter()
    indexes = {source: PositionIndex(ivs) for source, ivs in positions.items()}
    canonical = {source: {canon(v) for v in vids} for source, vids in excluded_vids.items()}
    for vid, c in candidates.items():
        reasons = [s for s, idx in sorted(indexes.items()) if idx.overlaps(*candidate_span(c))]
        reasons += [s for s, vids in sorted(canonical.items()) if canon(vid) in vids]
        if reasons:
            c.development = "+".join(reasons)
            for r in reasons:
                counts[r] += 1
    return counts


def exposure_class(c: Candidate) -> dict:
    """Novelty and the re-entry pattern for one candidate."""
    base_others = c.seen & {"E-dev", "T0"}
    return {
        "novel": not c.seen and c.development is None,
        # The re-entry pattern is defined on the three base sources only (E-dev, T0, E-hist), so
        # that adding archives cannot change whether they are required; every other set-A rule applies
        # (development material; the caller tallies only identity-resolved candidates).
        "novel_without_e_hist": not base_others and c.development is None,
        "re_entry": "E-hist" in c.seen and not base_others and c.development is None,
        "absent_all": not c.present,
        "upgrade": bool(c.seen) and isinstance(c.record.stars, int) and c.max_exposure_stars < c.record.stars,
    }


__all__ = ["BASE_EXPOSURE", "Candidate", "ExposureIndex", "InputFormatError", "PositionIndex", "T1Universe",
           "assemble_t1", "candidate_span", "check_info", "exposure_class", "key_span", "mark_development",
           "read_positions", "read_replacements", "read_submitters", "read_variant_summary", "read_vcf"]
