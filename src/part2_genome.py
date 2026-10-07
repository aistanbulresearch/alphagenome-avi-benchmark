"""Genome access and allele handling for the Part 2 harmonisation: a .fai-indexed FASTA reader, reverse
complement, left-normalisation and variant types."""
from __future__ import annotations

from pathlib import Path


_COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

def build_fai(fasta: Path, fai: Path) -> None:
    """Write a samtools-style .fai (name, length, offset, bases per line, bytes per line)."""
    out, name, length, offset, lb, lw = [], None, 0, 0, 0, 0
    pos = 0
    with open(fasta, "rb") as h:
        for line in h:
            if line.startswith(b">"):
                if name is not None:
                    out.append(f"{name}\t{length}\t{offset}\t{lb}\t{lw}")
                name, length, lb, lw = line[1:].split()[0].decode(), 0, 0, 0
                offset = pos + len(line)
            else:
                seq = line.rstrip(b"\r\n")
                if lb == 0:
                    lb, lw = len(seq), len(line)
                length += len(seq)
            pos += len(line)
    if name is not None:
        out.append(f"{name}\t{length}\t{offset}\t{lb}\t{lw}")
    fai.write_text("\n".join(out) + "\n", encoding="ascii")


class Fasta:
    """Random access to an uncompressed FASTA through its .fai (built next to it when absent)."""

    def __init__(self, path: Path | str, fai: Path | str | None = None):
        self.path = Path(path)
        fai = Path(fai) if fai else Path(str(self.path) + ".fai")
        if not fai.exists():
            build_fai(self.path, fai)
        self.index = {}
        for line in fai.read_text(encoding="ascii").splitlines():
            name, length, offset, lb, lw = line.split("\t")[:5]
            self.index[name] = (int(length), int(offset), int(lb), int(lw))
        self.handle = open(self.path, "rb")

    def _name(self, chrom: str) -> str:
        c = str(chrom)
        for cand in (c, "chr" + c, c[3:] if c.startswith("chr") else c):
            if cand in self.index:
                return cand
        raise KeyError(chrom)

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """Bases start..end (1-based, inclusive), upper case; '' outside the sequence."""
        try:
            length, offset, lb, lw = self.index[self._name(chrom)]
        except KeyError:
            return ""
        if start < 1 or end > length or end < start:
            return ""
        first = offset + (start - 1) // lb * lw + (start - 1) % lb
        last = offset + (end - 1) // lb * lw + (end - 1) % lb
        self.handle.seek(first)
        raw = self.handle.read(last - first + 1)
        return raw.replace(b"\n", b"").replace(b"\r", b"").decode("ascii").upper()


def revcomp(seq: str) -> str:
    return "".join(_COMP[b] for b in reversed(seq.upper()))


def normalize(fetch, chrom: str, pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Left-align and trim an allele (VCF convention; one anchor base kept for indels)."""
    ref, alt = ref.upper(), alt.upper()
    while True:
        if ref and alt and ref[-1] == alt[-1] and (len(ref) > 1 or len(alt) > 1):
            ref, alt = ref[:-1], alt[:-1]
        if not ref or not alt:
            base = fetch(chrom, pos - 1, pos - 1)
            if not base:
                break
            ref, alt, pos = base + ref, base + alt, pos - 1
            continue
        if ref[-1] != alt[-1] or (len(ref) == 1 and len(alt) == 1):
            break
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def variant_type(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(alt) == 1 and ref[0] == alt[0]:
        return "deletion"
    if len(ref) == 1 and alt[0] == ref[0]:
        return "insertion"
    return "delins"


def _chrom(c) -> str:
    c = str(c).strip()
    c = c[:-2] if c.endswith(".0") else c
    return c[3:] if c.lower().startswith("chr") else c


def _row(dataset, vid, chrom, pos, ref, alt, effect, measure, label, rule, gene, note, subset, use) -> dict:
    return {"dataset": dataset, "variant_id_in_source": str(vid), "chrom": _chrom(chrom), "pos": int(pos),
            "ref": str(ref).upper(), "alt": str(alt).upper(), "effect_size": effect, "effect_measure": measure,
            "splice_altering": label, "threshold_rule": rule, "gene": gene, "region_note": note,
            "subset": subset, "use": use}


def _label(flag: bool, intermediate: bool = False) -> str:
    return "yes" if flag else ("intermediate" if intermediate else "no")


def _true(v) -> bool:
    return v is True or (isinstance(v, str) and v.strip().upper() in ("TRUE", "YES")) or v == 1
