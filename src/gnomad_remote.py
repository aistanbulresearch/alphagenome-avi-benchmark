"""Read gnomAD v4.1 allele frequencies for a list of variants without downloading the full VCFs.

The joint (exomes + genomes) sites VCFs are bgzip-compressed with a tabix index, and the public bucket
serves byte ranges. For each variant the tabix index gives the compressed blocks that can hold its
position; only those blocks are fetched and decompressed. Standard library only.

Fields kept per variant (joint file): FILTER, AC_joint, AN_joint, AF_joint, fafmax_faf95_max_joint,
and whether the variant is present in the exome and genome call sets.
"""
from __future__ import annotations

import gzip
import struct
import urllib.request
import zlib
from dataclasses import dataclass
from typing import Callable, Iterable, Optional

BASE = "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint"
FIELDS = ("AC_joint", "AN_joint", "AF_joint", "fafmax_faf95_max_joint", "AC_exomes", "AC_genomes")
LINEAR_SHIFT = 14
MAX_BLOCK = 65536


def vcf_url(chrom: str) -> str:
    return f"{BASE}/gnomad.joint.v4.1.sites.chr{chrom}.vcf.bgz"


def fetch_range(url: str, start: int, end: int, timeout: int = 120) -> bytes:
    """Bytes start..end-1 of url (HTTP Range)."""
    req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end - 1}",
                                               "User-Agent": "alphagenome-avi-benchmark/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


@dataclass
class TabixIndex:
    names: list
    bins: list          # per reference: {bin: [(beg_voffset, end_voffset), ...]}
    linear: list        # per reference: [voffset, ...]

    def ref(self, name: str) -> int:
        return self.names.index(name)


def parse_tbi(data: bytes) -> TabixIndex:
    """Parse a tabix index (the .tbi bytes, still bgzip-compressed)."""
    raw = gzip.decompress(data)
    if raw[:4] != b"TBI\x01":
        raise ValueError("not a tabix index")
    n_ref = struct.unpack_from("<i", raw, 4)[0]
    l_nm = struct.unpack_from("<i", raw, 32)[0]
    names = [n.decode() for n in raw[36:36 + l_nm].split(b"\x00") if n]
    off = 36 + l_nm
    bins, linear = [], []
    for _ in range(n_ref):
        n_bin = struct.unpack_from("<i", raw, off)[0]
        off += 4
        ref_bins = {}
        for _ in range(n_bin):
            b, n_chunk = struct.unpack_from("<Ii", raw, off)
            off += 8
            chunks = list(struct.iter_unpack("<QQ", raw[off:off + 16 * n_chunk]))
            off += 16 * n_chunk
            ref_bins[b] = chunks
        n_intv = struct.unpack_from("<i", raw, off)[0]
        off += 4
        ioff = list(struct.unpack_from(f"<{n_intv}Q", raw, off))
        off += 8 * n_intv
        bins.append(ref_bins)
        linear.append(ioff)
    return TabixIndex(names, bins, linear)


def reg2bins(beg: int, end: int) -> list:
    """Bins overlapping [beg, end) (0-based), as in the SAM specification."""
    end -= 1
    out = [0]
    for shift, offset in ((26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)):
        out.extend(range(offset + (beg >> shift), offset + (end >> shift) + 1))
    return out


def chunks_for(index: TabixIndex, chrom_name: str, pos: int) -> list:
    """Merged (beg_voffset, end_voffset) chunks that can hold 1-based position `pos`."""
    rid = index.ref(chrom_name)
    beg = pos - 1
    lin = index.linear[rid]
    min_off = lin[min(beg >> LINEAR_SHIFT, len(lin) - 1)] if lin else 0
    chunks = sorted(c for b in reg2bins(beg, beg + 1) for c in index.bins[rid].get(b, ()) if c[1] > min_off)
    merged = []
    for c in chunks:
        if merged and c[0] <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], c[1]))
        else:
            merged.append(c)
    return merged


def inflate_blocks(data: bytes) -> list:
    """Decompress consecutive BGZF blocks; returns [(compressed_offset, text_bytes), ...]."""
    out, off = [], 0
    while off + 18 <= len(data):
        if data[off:off + 4] != b"\x1f\x8b\x08\x04":
            raise ValueError("not a BGZF block")
        bsize = struct.unpack_from("<H", data, off + 16)[0] + 1
        if off + bsize > len(data):
            break
        out.append((off, zlib.decompress(data[off + 18:off + bsize - 8], -15)))
        off += bsize
    return out


def read_chunk(url: str, chunk: tuple, fetch: Callable = fetch_range) -> bytes:
    """Uncompressed text of one tabix chunk."""
    cbeg, cend = chunk[0] >> 16, chunk[1] >> 16
    ubeg, uend = chunk[0] & 0xFFFF, chunk[1] & 0xFFFF
    data = fetch(url, cbeg, cend + MAX_BLOCK + 18)
    text = bytearray()
    for off, block in inflate_blocks(data):
        if off > cend - cbeg:
            break
        lo = ubeg if off == 0 else 0
        hi = uend if off == cend - cbeg else len(block)
        text += block[lo:hi]
    return bytes(text)


def parse_info(info: str) -> dict:
    out = {}
    for item in info.split(";"):
        key, _, value = item.partition("=")
        if key in FIELDS:
            out[key] = value
    return out


def lookup(url: str, index: TabixIndex, chrom: str, pos: int, alleles: Iterable[tuple],
           fetch: Callable = fetch_range) -> dict:
    """{(ref, alt): record or None} for alleles at one position; record holds FILTER and FIELDS."""
    wanted = {a: None for a in alleles}
    name = f"chr{chrom}"
    for chunk in chunks_for(index, name, pos):
        for line in read_chunk(url, chunk, fetch).decode("utf-8", "replace").split("\n"):
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t", 8)
            if len(parts) < 8 or parts[0] != name:
                continue
            p = int(parts[1])
            if p > pos:
                break
            if p == pos and (parts[3], parts[4]) in wanted:
                rec = parse_info(parts[7])
                rec["FILTER"] = parts[6]
                wanted[(parts[3], parts[4])] = rec
    return wanted


BGZF_MAGIC = b"\x1f\x8b\x08\x04"
PROBE = 24_576            # compressed bytes read to find the first complete record after an offset
SCAN = 98_304             # compressed bytes read per step of the final scan


def block_starts(data: bytes) -> list:
    """Offsets in `data` where a BGZF block starts (magic, 'BC' subfield and a consistent size)."""
    out, i = [], data.find(BGZF_MAGIC)
    while i != -1 and i + 18 <= len(data):
        if data[i + 12:i + 16] == b"BC\x02\x00":
            out.append(i)
            i = data.find(BGZF_MAGIC, i + struct.unpack_from("<H", data, i + 16)[0] + 1)
        else:
            i = data.find(BGZF_MAGIC, i + 1)
    return out


def first_position(url: str, offset: int, name: str, fetch: Callable = fetch_range) -> Optional[int]:
    """Position of the first complete record that starts after the first block start at or after
    `offset`, or None if none is found in PROBE bytes."""
    data = fetch(url, offset, offset + PROBE)
    starts = block_starts(data)
    if not starts:
        return None
    text = b"".join(block for _, block in inflate_blocks(data[starts[0]:]))
    for line in text.split(b"\n")[1:-1]:
        parts = line.split(b"\t", 2)
        if len(parts) > 2 and parts[0].decode() == name:
            return int(parts[1])
    return None


def scan_from(url: str, voffset: int, name: str, pos: int, fetch: Callable = fetch_range,
              skip_partial: bool = True) -> list:
    """Records at `pos`, reading forward from `voffset` until a record beyond `pos`. With `skip_partial`
    the text before the first newline is dropped (an offset that is not a record start)."""
    off, within = voffset >> 16, voffset & 0xFFFF
    carry, found, first = b"", [], True
    while True:
        data = fetch(url, off, off + SCAN)
        starts = block_starts(data)
        if not starts:
            return found
        blocks = inflate_blocks(data[starts[0]:])
        if not blocks:
            return found
        text = b"".join(b for _, b in blocks)
        if first:
            text = text[within:] if starts[0] == 0 else text
            if skip_partial:
                text = text.split(b"\n", 1)[1] if b"\n" in text else b""
            first = False
        text = carry + text
        lines = text.split(b"\n")
        carry = lines.pop()
        for line in lines:
            parts = line.decode("utf-8", "replace").split("\t", 8)
            if len(parts) < 2 or parts[0] != name:          # VCF has 8+ columns, CADD has 6
                continue
            p = int(parts[1])
            if p > pos:
                return found
            if p == pos:
                found.append(parts)
        last_off = starts[0] + blocks[-1][0]
        off = off + last_off + struct.unpack_from("<H", data, last_off + 16)[0] + 1


def search(url: str, index: TabixIndex, chrom: str, pos: int, fetch: Callable = fetch_range,
           max_probes: int = 30, contig: Optional[str] = None) -> list:
    """Records at 1-based `pos` of any tab-separated, tabix-indexed file with the sequence name in the
    first column and the position in the second (VCF, CADD TSV). `contig` is the sequence name as the
    file writes it (default 'chr' + chrom). Within the tabix chunk for the position, interpolation and
    bisection over compressed offsets find a block whose first complete record lies before `pos`; the
    records are then read forward from there."""
    name = contig or f"chr{chrom}"
    found = []
    for chunk in chunks_for(index, name, pos):
        lo, hi = chunk[0] >> 16, chunk[1] >> 16
        start = None                                   # compressed offset with first record < pos
        p_lo, p_hi = None, None
        for _ in range(max_probes):
            if hi - lo <= SCAN:
                break
            if p_lo is not None and p_hi is not None and p_hi > p_lo:
                guess = lo + int((pos - p_lo) / (p_hi - p_lo) * (hi - lo)) - PROBE // 2
            else:
                guess = (lo + hi) // 2
            margin = (hi - lo) // 8                        # keep each step a real cut of the interval
            guess = max(lo + margin, min(guess, hi - margin))
            p = first_position(url, guess, name, fetch)
            if p is None or p >= pos:
                hi, p_hi = guess, p if p is not None else p_hi
            else:
                lo, p_lo, start = guess, p, guess
        if start is None:
            found += scan_from(url, chunk[0], name, pos, fetch, skip_partial=False)
        else:
            found += scan_from(url, start << 16, name, pos, fetch, skip_partial=True)
    return found


def lookup_fast(url: str, index: TabixIndex, chrom: str, pos: int, alleles: Iterable[tuple],
                fetch: Callable = fetch_range) -> dict:
    """As `lookup`, reading only the blocks near the position."""
    wanted = {a: None for a in alleles}
    for parts in search(url, index, chrom, pos, fetch):
        if (parts[3], parts[4]) in wanted:
            rec = parse_info(parts[7])
            rec["FILTER"] = parts[6]
            wanted[(parts[3], parts[4])] = rec
    return wanted


def as_float(value: Optional[str]) -> Optional[float]:
    try:
        return float(value) if value not in (None, "", ".") else None
    except ValueError:
        return None


__all__ = ["vcf_url", "fetch_range", "parse_tbi", "reg2bins", "chunks_for", "inflate_blocks", "read_chunk",
           "lookup", "parse_info", "as_float", "TabixIndex"]
