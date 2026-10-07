"""Per-base values from a remote bigWig file (UCSC format) by HTTP byte ranges.

Reads the header, the chromosome B+ tree and the R-tree index once, then for each requested position
the data block that covers it (zlib-compressed; bedGraph, variableStep or fixedStep sections).
Standard library only. Format: Kent et al. 2010, Bioinformatics 26:2204, and the UCSC source.
"""
from __future__ import annotations

import struct
import zlib
from bisect import bisect_right
from typing import Callable, Optional

from gnomad_remote import fetch_range

BIGWIG_MAGIC = 0x888FFC26
CHROM_TREE_MAGIC = 0x78CA8C91
RTREE_MAGIC = 0x2468ACE0


class BigWig:
    def __init__(self, url: str, fetch: Callable = fetch_range):
        self.url, self.fetch = url, fetch
        head = fetch(url, 0, 64)
        (magic, _version, _zooms, chrom_tree, _data, index, _fields, _defined, _auto, _total,
         uncompress) = struct.unpack_from("<IHHQQQHHQQI", head, 0)
        if magic != BIGWIG_MAGIC:
            raise ValueError("not a bigWig file")
        self.compressed = uncompress > 0
        self.chroms = self._read_chrom_tree(chrom_tree)
        self.index_offset = index
        self._node_cache: dict = {}

    def _read_chrom_tree(self, offset: int) -> dict:
        h = self.fetch(self.url, offset, offset + 32)
        magic, block_size, key_size, val_size, count = struct.unpack_from("<IIIIQ", h, 0)
        if magic != CHROM_TREE_MAGIC:
            raise ValueError("bad chromosome tree")
        out = {}

        def node(off: int) -> None:
            raw = self.fetch(self.url, off, off + 4 + block_size * (key_size + 8))
            is_leaf, _, n = struct.unpack_from("<BBH", raw, 0)
            p = 4
            for _ in range(n):
                key = raw[p:p + key_size].rstrip(b"\x00").decode()
                if is_leaf:
                    cid, size = struct.unpack_from("<II", raw, p + key_size)
                    out[key] = (cid, size)
                    p += key_size + 8
                else:
                    child = struct.unpack_from("<Q", raw, p + key_size)[0]
                    node(child)
                    p += key_size + 8
        node(offset + 32)
        return out

    def _node(self, off: int) -> tuple:
        if off not in self._node_cache:
            raw = self.fetch(self.url, off, off + 4)
            is_leaf, _, n = struct.unpack_from("<BBH", raw, 0)
            size = 32 if is_leaf else 24
            body = self.fetch(self.url, off + 4, off + 4 + n * size)
            items = []
            for i in range(n):
                sc, sb, ec, eb = struct.unpack_from("<IIII", body, i * size)
                if is_leaf:
                    data_off, data_size = struct.unpack_from("<QQ", body, i * size + 16)
                    items.append((sc, sb, ec, eb, data_off, data_size))
                else:
                    items.append((sc, sb, ec, eb, struct.unpack_from("<Q", body, i * size + 16)[0], None))
            self._node_cache[off] = (is_leaf, items)
        return self._node_cache[off]

    def blocks(self, cid: int, start: int, end: int) -> list:
        """Data blocks (offset, size) overlapping [start, end) on chromosome id cid (0-based)."""
        out = []

        def overlaps(it):
            return (it[0], it[1]) < (cid, end) and (it[2], it[3]) > (cid, start)

        def walk(off):
            is_leaf, items = self._node(off)
            for it in items:
                if overlaps(it):
                    if is_leaf:
                        out.append((it[4], it[5]))
                    else:
                        walk(it[4])
        walk(self.index_offset + 48)
        return out

    def block_values(self, off: int, size: int) -> list:
        """[(start, end, value)] 0-based half-open intervals of one data block."""
        raw = self.fetch(self.url, off, off + size)
        data = zlib.decompress(raw) if self.compressed else raw
        out, p = [], 0
        while p + 24 <= len(data):
            _cid, start, _end, step, span, typ, _r, n = struct.unpack_from("<IIIIIBBH", data, p)
            p += 24
            for i in range(n):
                if typ == 1:                                   # bedGraph
                    s, e, v = struct.unpack_from("<IIf", data, p)
                    p += 12
                elif typ == 2:                                 # variableStep
                    s, v = struct.unpack_from("<If", data, p)
                    e = s + span
                    p += 8
                else:                                          # fixedStep
                    v = struct.unpack_from("<f", data, p)[0]
                    s = start + i * step
                    e = s + span
                    p += 4
                out.append((s, e, v))
        return out

    def value(self, chrom: str, pos: int, cache: Optional[dict] = None) -> Optional[float]:
        """Value at 1-based position, or None if not covered."""
        if chrom not in self.chroms:
            return None
        cid = self.chroms[chrom][0]
        for off, size in self.blocks(cid, pos - 1, pos):
            key = (off, size)
            vals = cache.get(key) if cache is not None else None
            if vals is None:
                vals = self.block_values(off, size)
                if cache is not None:
                    cache[key] = vals
            starts = [s for s, _, _ in vals]
            i = bisect_right(starts, pos - 1) - 1
            if i >= 0 and vals[i][0] <= pos - 1 < vals[i][1]:
                return float(vals[i][2])
        return None


__all__ = ["BigWig"]
