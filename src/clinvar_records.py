"""Record-level rules for ClinVar variant_summary rows and VCF lines.

Every function is pure and is tested on synthetic rows. Field names are taken from each file's header.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Iterable, Optional, Union

# eligibility at T1 is an exact match of the whole field.
STRICT_PLP = ("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic")
STRICT_BLB = ("Benign", "Likely benign", "Benign/Likely benign")
# Exposure tokens after splitting on ';', ',', '|' and '/', compared without regard to case.
SEEN_TOKENS = frozenset({"pathogenic", "likely pathogenic", "benign", "likely benign"})
# Every other classification term ClinVar uses for germline records (lower case, after splitting). A
# token outside both sets is unknown and makes the record count as seen (the conservative direction).
KNOWN_OTHER_TOKENS = frozenset({
    "uncertain significance", "conflicting classifications of pathogenicity",
    "conflicting interpretations of pathogenicity", "not provided", "drug response", "risk factor",
    "association", "association not found", "protective", "affects", "other", "low penetrance",
    "established risk allele", "likely risk allele", "uncertain risk allele", "confers sensitivity",
    "no classification for the single variant", "no classifications from unflagged records",
    "no interpretation for the single variant", "-", "."})   # '-' and '.' mean no classification
_SPLIT = re.compile(r"[;,|/]")
_INCL_START = re.compile(r"^(\d+):(.*)$")


def incl_pairs(value: str) -> list:
    """(VariationID, classification) pairs of a CLNSIGINCL value. ClinVar uses '|' both between pairs
    and between the terms of one classification: a token that begins '<digits>:' opens a pair, any other
    token without a colon continues the classification of the previous pair. A value whose first token
    is not a pair, a colon-bearing token that is not a pair, an empty token or classification, or
    whitespace inside a token is malformed (ValueError)."""
    pairs: list = []
    for token in value.split("|"):
        if not token or any(ch.isspace() for ch in token):
            raise ValueError(f"malformed CLNSIGINCL token {token!r}")
        m = _INCL_START.match(token)
        if m:
            if not m.group(2):
                raise ValueError(f"CLNSIGINCL pair without a classification: {token!r}")
            pairs.append([m.group(1), m.group(2)])
        elif ":" in token:                                  # a colon only ever opens a pair
            raise ValueError(f"malformed CLNSIGINCL pair {token!r}")
        elif pairs:
            pairs[-1][1] += "|" + token
        else:
            raise ValueError(f"CLNSIGINCL value does not begin with a pair: {token!r}")
    return [tuple(p) for p in pairs]

STAR_MAP = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, single submitter": 1,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "no assertion criteria provided": 0,
    "no classification provided": 0,
    "no classification for the single variant": 0,
    "no interpretation for the single variant": 0,
}
UNKNOWN = "unknown"
Stars = Union[int, str]

# chromosome groups from the AVI preprint (Methods, Data Split).
CHROM_GROUPS = {
    "training": ("1", "4", "7", "8", "10", "13", "15"),
    "validation": ("2", "5", "11", "14", "17", "20", "22", "X"),
    "untouched": ("3", "6", "9", "12", "16", "18", "19", "21"),
}
_GROUP_OF = {c: g for g, chroms in CHROM_GROUPS.items() for c in chroms}
SCOPES = {"S1": ("untouched",), "S2": ("untouched", "validation"),
          "S3": ("untouched", "validation", "training")}
_BASES = re.compile(r"^[ACGTN]+$")
MAX_INDEL = 10
REQUIRED_VS_FIELDS = ("#AlleleID", "Name", "GeneID", "ClinicalSignificance", "ReviewStatus", "Assembly",
                      "Chromosome", "Start", "Stop", "VariationID", "PositionVCF", "ReferenceAlleleVCF",
                      "AlternateAlleleVCF")


def strict_group(clinical_significance: Optional[str]) -> Optional[str]:
    """'P/LP' or 'B/LB' if the whole field is one of the strict values, else None."""
    if clinical_significance in STRICT_PLP:
        return "P/LP"
    if clinical_significance in STRICT_BLB:
        return "B/LB"
    return None


def _tokens(text: Optional[str]) -> Iterable[str]:
    for token in _SPLIT.split(text or ""):
        yield token.replace("_", " ").strip().lower()


def unknown_tokens(clinical_significance: Optional[str]) -> list[str]:
    """Non-empty tokens that are neither strict nor a known other ClinVar classification term."""
    return [t for t in _tokens(clinical_significance) if t and t not in SEEN_TOKENS and t not in KNOWN_OTHER_TOKENS]


def strict_present(clinical_significance: Optional[str]) -> bool:
    """True if a strict term (P, LP, B, LB) is among the tokens."""
    return any(t in SEEN_TOKENS for t in _tokens(clinical_significance))


def seen_strict(clinical_significance: Optional[str]) -> bool:
    """Exposure rule: any token equals a P, LP, B or LB classification, or any token is unknown."""
    tokens = [t for t in _tokens(clinical_significance) if t]
    return any(t in SEEN_TOKENS or t not in KNOWN_OTHER_TOKENS for t in tokens)


def vcf_seen_strict(clnsig: Optional[str], clnsigincl: Optional[str]) -> bool:
    """Exposure rule for the VCF: CLNSIG, or any classification in CLNSIGINCL ('VariationID:class')."""
    if seen_strict(clnsig):
        return True
    if not clnsigincl:
        return False
    return any(seen_strict(classification) for _, classification in incl_pairs(clnsigincl))


def stars(review_status: Optional[str]) -> Stars:
    """Stars from a variant_summary ReviewStatus or a VCF CLNREVSTAT value; 'unknown' otherwise."""
    if review_status is None:
        return UNKNOWN
    key = review_status.replace("_", " ").strip()
    key = re.sub(r"\s*,\s*", ", ", key)
    return STAR_MAP.get(key, UNKNOWN)


def chromosome_group(chrom: Optional[str]) -> str:
    """training, validation or untouched; 'excluded' for Y, MT and anything else."""
    return _GROUP_OF.get(canonical_chrom(chrom), "excluded")


def in_scope(group: str, scope: str) -> bool:
    return group in SCOPES[scope]


def canonical_chrom(chrom: Optional[str]) -> str:
    """Upper case, without a 'chr' prefix ('chrx' -> 'X')."""
    c = (chrom or "").strip().upper()
    return c[3:] if c.startswith("CHR") else c


def normalize_allele(chrom: str, pos: Union[int, str], ref: str, alt: str) -> Optional[str]:
    """Key 'chrom:pos:ref:alt' after trimming shared trailing then leading bases (one base kept).

    Left alignment needs the reference sequence and is not done; both sources are ClinVar's own
    VCF-style representations, which ClinVar documents as left-shifted. The chromosome is upper-cased
    with any 'chr' prefix removed. Returns None when the allele cannot be represented.
    """
    try:
        p = int(pos)
    except (TypeError, ValueError):
        return None
    if (ref or "").strip().lower() in ("", "na", "-") or (alt or "").strip().lower() in ("", "na", "-"):
        return None
    ref, alt = ref.strip().upper(), alt.strip().upper()
    c = canonical_chrom(chrom)
    if p < 1 or not c or not _BASES.match(ref) or not _BASES.match(alt) or ref == alt:
        return None
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, p = ref[1:], alt[1:], p + 1
    return f"{c}:{p}:{ref}:{alt}"


def variant_type(key: Optional[str]) -> str:
    """'snv'; 'small_indel' (pure insertion or deletion of 1-10 bp after a shared anchor base); 'other'."""
    if key is None:
        return "other"
    ref, alt = key.split(":")[2:4]
    if len(ref) == 1 and len(alt) == 1:
        return "snv"
    change = abs(len(ref) - len(alt))
    if 1 <= change <= MAX_INDEL and min(len(ref), len(alt)) == 1 and ref[0] == alt[0]:
        return "small_indel"
    return "other"


@dataclass(frozen=True)
class Record:
    """One allele-level ClinVar record from either file type."""

    source: str
    allele_id: Optional[str]
    variation_id: Optional[str]
    assembly: str
    chrom: str
    key: Optional[str]          # normalized GRCh38 allele, None otherwise
    name: str
    classification: str
    review_status: str
    stars: Stars
    seen: bool                  # strict for exposure
    gene: str
    key37: Optional[str] = None  # normalized GRCh37 allele (variant_summary GRCh37 rows only)
    strict_token: bool = False   # seen because of a strict term (not only an unknown term)
    span: Optional[tuple] = None  # (chrom, start, end) on `assembly`, kept even without an allele key
    unknown_terms: tuple = ()     # unique unknown classification terms (CLNSIG and every CLNSIGINCL classification)

    def identity(self) -> tuple:
        """Fields that must agree for duplicate rows of one VariationID to be collapsed."""
        return (self.allele_id, self.key, self.classification, self.review_status, self.name)


def placement(chrom: Optional[str], start, end=None, ref: Optional[str] = None) -> Optional[tuple]:
    """(chrom, start, end), 1-based inclusive; end from `end`, else from the reference length, else start."""
    c = canonical_chrom(chrom)
    try:
        s = int(start)
    except (TypeError, ValueError):
        return None
    if not c or s < 1:
        return None
    e = s
    try:
        e = max(s, int(end)) if end not in (None, "", "-", "na") else e
    except (TypeError, ValueError):
        pass
    if end in (None, "", "-", "na") and ref and _BASES.match(ref.strip().upper()):
        e = s + len(ref.strip()) - 1
    return c, s, e


def _clean_id(value: Optional[str]) -> Optional[str]:
    v = (value or "").strip()
    return v if v and v not in ("-1", "-", "na", ".") else None


def from_variant_summary(row: dict, source: str) -> Record:
    """Record from one variant_summary row given as {header field: value}."""
    assembly = row.get("Assembly", "").strip()
    chrom = row.get("Chromosome", "").strip()
    allele = normalize_allele(chrom, row.get("PositionVCF"), row.get("ReferenceAlleleVCF"),
                              row.get("AlternateAlleleVCF"))
    key = allele if assembly == "GRCh38" else None
    key37 = allele if assembly == "GRCh37" else None
    clinsig = row.get("ClinicalSignificance", "")
    review = row.get("ReviewStatus", "")
    gene = _clean_id(row.get("GeneID")) or "missing"
    span = placement(chrom, row.get("Start"), row.get("Stop")) or placement(
        chrom, row.get("PositionVCF"), ref=row.get("ReferenceAlleleVCF"))
    return Record(source, _clean_id(row.get("#AlleleID")), _clean_id(row.get("VariationID")), assembly,
                  chrom, key, row.get("Name", ""), clinsig, review, stars(review), seen_strict(clinsig), gene,
                  key37, strict_present(clinsig), span, tuple(dict.fromkeys(unknown_tokens(clinsig))))


def parse_info(info: str) -> dict:
    out = {}
    for item in info.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
        elif item:
            out[item] = ""
    return out


def from_vcf_line(line: str, source: str) -> list[Record]:
    """Records from one GRCh38 VCF data line: one per ALT allele; for ALT '.' one record without an
    allele key, so that its identifiers and classification still enter the exposure index."""
    fields = line.rstrip("\r\n").split("\t")
    chrom, pos, vid, ref, alts, info = fields[0], fields[1], fields[2], fields[3], fields[4], fields[7]
    tags = parse_info(info)
    clnsig, review = tags.get("CLNSIG", ""), tags.get("CLNREVSTAT", "")
    seen = vcf_seen_strict(clnsig, tags.get("CLNSIGINCL"))
    included = [c for _, c in incl_pairs(tags["CLNSIGINCL"])] if tags.get("CLNSIGINCL") else []
    strict = strict_present(clnsig) or any(strict_present(c) for c in included)
    unknown = tuple(dict.fromkeys(unknown_tokens(clnsig) + [t for c in included for t in unknown_tokens(c)]))
    gene = "missing"
    if tags.get("GENEINFO"):
        gene = "|".join(p.split(":", 1)[-1] for p in tags["GENEINFO"].split("|"))
    records = []
    for alt in alts.split(","):
        key = None if alt == "." else normalize_allele(chrom, pos, ref, alt)
        records.append(Record(source, _clean_id(tags.get("ALLELEID")), _clean_id(vid), "GRCh38", chrom,
                              key, tags.get("CLNHGVS", ""), clnsig,
                              review, stars(review), seen, gene, None, strict, placement(chrom, pos, ref=ref),
                              unknown))
    return records


__all__ = [
    "CHROM_GROUPS", "MAX_INDEL", "REQUIRED_VS_FIELDS", "Record", "canonical_chrom", "SCOPES", "STAR_MAP", "STRICT_BLB", "STRICT_PLP", "UNKNOWN",
    "KNOWN_OTHER_TOKENS", "chromosome_group", "from_variant_summary", "from_vcf_line", "in_scope",
    "normalize_allele", "placement", "unknown_tokens",
    "parse_info", "seen_strict", "stars", "strict_group", "variant_type", "vcf_seen_strict",
]
