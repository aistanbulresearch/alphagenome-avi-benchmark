"""GRCh37 -> GRCh38 interval conversion with a UCSC chain file.

Reads the chain file `data/derived/exposure/hg19ToHg38.over.chain.gz` (UCSC). A span is mapped
through every aligned block it overlaps, in every chain; the result is the list of GRCh38 spans of the
overlapping parts (several when the span is split or maps to several places). An empty list means the
span does not convert. Coordinates are 1-based inclusive at the interface; chain blocks are 0-based
half-open internally. Standard library only.
"""
from __future__ import annotations

import bisect
import gzip
from collections import defaultdict
from pathlib import Path

import clinvar_records as cr


class Lifter:
    def __init__(self, blocks: dict[str, list[tuple[int, int, str, int, str, int]]]):
        # per source chromosome: sorted (t_start, t_end, q_chrom, q_start, q_strand, q_size)
        self.blocks = {c: sorted(b) for c, b in blocks.items()}
        self.starts = {c: [b[0] for b in bl] for c, bl in self.blocks.items()}
        self.max_len = {c: max((b[1] - b[0] for b in bl), default=0) for c, bl in self.blocks.items()}

    def __call__(self, span: tuple[str, int, int]) -> list[tuple[str, int, int]]:
        chrom, start, end = cr.canonical_chrom(span[0]), span[1] - 1, span[2]     # to 0-based half-open
        if chrom not in self.blocks:
            return []
        out = []
        bl, starts = self.blocks[chrom], self.starts[chrom]
        i = max(bisect.bisect_left(starts, start - self.max_len[chrom]) - 1, 0)
        for t_start, t_end, q_chrom, q_start, strand, q_size in bl[i:]:
            if t_start >= end:
                break
            lo, hi = max(start, t_start), min(end, t_end)
            if lo >= hi:
                continue
            a, b = q_start + (lo - t_start), q_start + (hi - t_start)
            if strand == "-":
                a, b = q_size - b, q_size - a
            out.append((cr.canonical_chrom(q_chrom), a + 1, b))
        return out


def parse_chain(lines) -> Lifter:
    blocks = defaultdict(list)
    header = None
    for raw in lines:
        line = raw.strip()
        if not line:
            header = None
            continue
        parts = line.split()
        if parts[0] == "chain":
            # chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
            header = {"t": cr.canonical_chrom(parts[2]), "t_pos": int(parts[5]), "q": parts[7],
                      "q_size": int(parts[8]), "q_strand": parts[9], "q_pos": int(parts[10])}
            continue
        if header is None:
            raise ValueError("chain data line before a chain header")
        size = int(parts[0])
        blocks[header["t"]].append((header["t_pos"], header["t_pos"] + size, header["q"], header["q_pos"],
                                    header["q_strand"], header["q_size"]))
        if len(parts) == 3:
            header["t_pos"] += size + int(parts[1])
            header["q_pos"] += size + int(parts[2])
        else:
            header["t_pos"] += size
            header["q_pos"] += size
    return Lifter(blocks)


def load(path: Path) -> Lifter:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return parse_chain(handle)


__all__ = ["Lifter", "load", "parse_chain"]
