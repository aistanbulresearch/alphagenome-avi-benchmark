"""Records at chosen positions of a remote BGZF file sorted by position that has no tabix index.

A sparse map of (compressed offset, first position after that offset) is built once per file by probing
at a fixed byte step. A group of nearby target positions (a window) is then located by interpolation
between the two map entries around it, and read forward from there until a record beyond the window.
The file holds one sequence (contig); its name is in the first column and the position in the second.
"""
from __future__ import annotations

import bisect
import struct
from concurrent.futures import ThreadPoolExecutor
from typing import Callable, Optional

import gnomad_remote as gr

PROBE = 40_960        # compressed bytes per probe: more than two BGZF blocks of the SpliceAI2 files
SCAN = 98_304         # compressed bytes per step of a forward read


def first_pos(url: str, offset: int, name: str, fetch: Callable, probe: int = PROBE) -> Optional[int]:
    """Position of the first complete record after the first block start at or after `offset`."""
    data = fetch(url, offset, offset + probe)
    starts = gr.block_starts(data)
    if not starts:
        return None
    text = b"".join(block for _, block in gr.inflate_blocks(data[starts[0]:]))
    for line in text.split(b"\n")[1:-1]:
        parts = line.split(b"\t", 2)
        if len(parts) > 2 and parts[0].decode() == name:
            return int(parts[1])
    return None


def sparse_map(url: str, size: int, name: str, fetch: Callable, step: int, probe: int = PROBE,
               threads: int = 16) -> list:
    """[(offset, first position or None)] at offsets 0, step, 2*step, ... below `size`."""
    offsets = list(range(0, size, step))
    with ThreadPoolExecutor(threads) as pool:
        found = list(pool.map(lambda o: first_pos(url, o, name, fetch, probe), offsets))
    return list(zip(offsets, found))


def bracket(smap: list, size: int, pos: int) -> tuple:
    """(lo, p_lo, hi, p_hi): the last map offset whose first position is below `pos` (0 if none) and the
    first later offset whose first position is at or beyond it (the file size if none)."""
    known = [(o, p) for o, p in smap if p is not None]
    keys = [p for _, p in known]
    i = bisect.bisect_left(keys, pos)
    lo, p_lo = known[i - 1] if i > 0 else (0, None)
    hi, p_hi = known[i] if i < len(known) else (size, None)
    return lo, p_lo, hi, p_hi


def locate(url: str, name: str, pos: int, lo: int, p_lo: Optional[int], hi: int, p_hi: Optional[int],
           fetch: Callable, probe: int = PROBE, scan: int = SCAN, max_probes: int = 30) -> int:
    """A compressed offset whose first record lies before `pos` (or 0), found by interpolation."""
    for _ in range(max_probes):
        if hi - lo <= scan:
            break
        if p_lo is not None and p_hi is not None and p_hi > p_lo:
            guess = lo + int((pos - p_lo) / (p_hi - p_lo) * (hi - lo)) - probe // 2
        else:
            guess = (lo + hi) // 2
        margin = (hi - lo) // 8
        guess = max(lo + margin, min(guess, hi - margin))
        p = first_pos(url, guess, name, fetch, probe)
        if p is None or p >= pos:
            hi, p_hi = guess, p if p is not None else p_hi
        else:
            lo, p_lo = guess, p
    return lo


def scan_window(url: str, start: int, name: str, wanted: set, end: int, fetch: Callable,
                scan: int = SCAN) -> dict:
    """{position: [record fields]} for records at `wanted` positions, read forward from the first block
    start at or after `start` (its first, possibly partial, line skipped) until a record beyond `end`."""
    off, carry, first, out = start, b"", True, {}
    while True:
        data = fetch(url, off, off + scan)
        starts = gr.block_starts(data)
        if not starts:
            return out
        blocks = gr.inflate_blocks(data[starts[0]:])
        if not blocks:
            return out
        text = b"".join(b for _, b in blocks)
        if first:
            text = text.split(b"\n", 1)[1] if b"\n" in text else b""
            first = False
        lines = (carry + text).split(b"\n")
        carry = lines.pop()
        for line in lines:
            parts = line.decode("utf-8", "replace").split("\t")
            if len(parts) < 2 or parts[0] != name:
                continue
            p = int(parts[1])
            if p > end:
                return out
            if p in wanted:
                out.setdefault(p, []).append(parts)
        last = starts[0] + blocks[-1][0]
        off = off + last + struct.unpack_from("<H", data, last + 16)[0] + 1


def windows(positions: list, gap: int) -> list:
    """Sorted distinct positions merged into (first, last) windows wherever neighbours are within `gap`."""
    out = []
    for p in sorted(set(positions)):
        if out and p - out[-1][1] <= gap:
            out[-1] = (out[-1][0], p)
        else:
            out.append((p, p))
    return out


def fetch_positions(url: str, size: int, smap: list, name: str, positions: list, fetch: Callable,
                    probe: int = PROBE, scan: int = SCAN, gap: int = 200, threads: int = 1) -> dict:
    """{position: [record fields]} for every target position found in the file; windows run in parallel."""
    wanted = set(positions)

    def one(window):
        a, b = window
        lo, p_lo, hi, p_hi = bracket(smap, size, a)
        start = locate(url, name, a, lo, p_lo, hi, p_hi, fetch, probe, scan)
        own = {p for p in wanted if a <= p <= b}      # each window keeps only its own positions
        return scan_window(url, start, name, own, b, fetch, scan)

    out = {}
    with ThreadPoolExecutor(threads) as pool:
        for found in pool.map(one, windows(positions, gap)):
            out.update(found)
    return out


__all__ = ["first_pos", "sparse_map", "bracket", "locate", "scan_window", "windows", "fetch_positions"]
