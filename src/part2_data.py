"""Part 2: experimental splicing measurements, downloaded and harmonised to GRCh38.

Placeholder header; replaced below.
"""
from __future__ import annotations

import io
import math
import re
import struct
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

# ---------------------------------------------------------------------------------------------------
# Spreadsheet readers (standard library only): .xlsx (Office Open XML) and .xls (BIFF8 in an OLE2
# compound file). Each returns {sheet name: list of rows}, a row being a list of cell values
# (str, float, bool or None).
# ---------------------------------------------------------------------------------------------------

_NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
       "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
       "p": "http://schemas.openxmlformats.org/package/2006/relationships"}


def _col_index(ref: str) -> int:
    n = 0
    for ch in ref:
        if ch.isalpha():
            n = n * 26 + (ord(ch.upper()) - 64)
        else:
            break
    return n - 1


def read_xlsx(path: Path | str | bytes) -> dict[str, list[list]]:
    data = path if isinstance(path, bytes) else Path(path).read_bytes()
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        names = set(zf.namelist())
        shared: list[str] = []
        if "xl/sharedStrings.xml" in names:
            root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            for si in root.findall("m:si", _NS):
                shared.append("".join(t.text or "" for t in si.iter(f"{{{_NS['m']}}}t")))
        rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        target = {r.get("Id"): r.get("Target") for r in rels.findall("p:Relationship", _NS)}
        wb = ET.fromstring(zf.read("xl/workbook.xml"))
        out: dict[str, list[list]] = {}
        for sheet in wb.find("m:sheets", _NS).findall("m:sheet", _NS):
            rid = sheet.get(f"{{{_NS['r']}}}id")
            t = target[rid].lstrip("/")
            member = t if t.startswith("xl/") else "xl/" + t
            root = ET.fromstring(zf.read(member))
            rows: dict[int, dict[int, object]] = {}
            for row in root.iter(f"{{{_NS['m']}}}row"):
                r_idx = int(row.get("r")) - 1
                cells = rows.setdefault(r_idx, {})
                for c in row.findall("m:c", _NS):
                    ref, typ = c.get("r"), c.get("t")
                    v = c.find("m:v", _NS)
                    if typ == "inlineStr":
                        val = "".join(t.text or "" for t in c.iter(f"{{{_NS['m']}}}t"))
                    elif v is None or v.text is None:
                        continue
                    elif typ == "s":
                        val = shared[int(v.text)]
                    elif typ in ("str", "e"):
                        val = v.text
                    elif typ == "b":
                        val = v.text == "1"
                    else:
                        val = float(v.text)
                    cells[_col_index(ref)] = val
            out[sheet.get("name")] = _dense(rows)
        return out


def _dense(rows: dict[int, dict[int, object]]) -> list[list]:
    if not rows:
        return []
    n_rows = max(rows) + 1
    n_cols = max((max(c) + 1 for c in rows.values() if c), default=0)
    grid = [[None] * n_cols for _ in range(n_rows)]
    for r, cells in rows.items():
        for c, v in cells.items():
            grid[r][c] = v
    return grid


def _cfb_stream(data: bytes, wanted: tuple[str, ...]) -> bytes:
    """Return the first stream named in `wanted` from an OLE2 compound file."""
    if data[:8] != bytes.fromhex("D0CF11E0A1B11AE1"):
        raise ValueError("not an OLE2 compound file")
    sec_size = 1 << struct.unpack_from("<H", data, 30)[0]
    mini_size = 1 << struct.unpack_from("<H", data, 32)[0]
    n_fat, dir_start = struct.unpack_from("<II", data, 44)
    mini_cutoff, minifat_start, n_minifat, difat_start, n_difat = struct.unpack_from("<IIIII", data, 56)

    def sector(i: int) -> bytes:
        off = (i + 1) * sec_size
        return data[off:off + sec_size]

    difat = list(struct.unpack_from("<109I", data, 76))
    nxt = difat_start
    for _ in range(n_difat):
        block = struct.unpack(f"<{sec_size // 4}I", sector(nxt))
        difat.extend(block[:-1])
        nxt = block[-1]
    fat: list[int] = []
    for s in difat[:n_fat]:
        fat.extend(struct.unpack(f"<{sec_size // 4}I", sector(s)))

    def chain(start: int, table: list[int]) -> list[int]:
        out, s = [], start
        while s < 0xFFFFFFFA:
            out.append(s)
            s = table[s]
        return out

    directory = b"".join(sector(s) for s in chain(dir_start, fat))
    entries = []
    for off in range(0, len(directory), 128):
        e = directory[off:off + 128]
        name_len = struct.unpack_from("<H", e, 64)[0]
        name = e[:max(name_len - 2, 0)].decode("utf-16-le", "replace")
        start, size = struct.unpack_from("<II", e, 116)
        entries.append((name, e[66], start, size))
    root_start, root_size = entries[0][2], entries[0][3]
    for name, typ, start, size in entries:
        if typ != 2 or name not in wanted:
            continue
        if size >= mini_cutoff:
            return b"".join(sector(s) for s in chain(start, fat))[:size]
        ministream = b"".join(sector(s) for s in chain(root_start, fat))[:root_size]
        minifat: list[int] = []
        for s in chain(minifat_start, fat):
            minifat.extend(struct.unpack(f"<{sec_size // 4}I", sector(s)))
        return b"".join(ministream[s * mini_size:(s + 1) * mini_size] for s in chain(start, minifat))[:size]
    raise ValueError(f"no stream named {wanted} in the compound file")


def _rk(value: int) -> float:
    if value & 2:
        num = float(value >> 2 if not value & 0x80000000 else (value >> 2) - (1 << 30))
    else:
        num = struct.unpack("<d", struct.pack("<Q", (value & 0xFFFFFFFC) << 32))[0]
    return num / 100 if value & 1 else num


class _Pieces:
    """Byte reader over a record body and its CONTINUE records (needed for the shared string table)."""

    def __init__(self, parts: list[bytes]):
        self.parts, self.i, self.pos = parts, 0, 0

    def _avail(self) -> int:
        return len(self.parts[self.i]) - self.pos

    def take(self, n: int) -> bytes:
        out = b""
        while n:
            if self._avail() == 0:
                self.i, self.pos = self.i + 1, 0
            k = min(n, self._avail())
            out += self.parts[self.i][self.pos:self.pos + k]
            self.pos += k
            n -= k
        return out

    def chars(self, count: int, wide: bool) -> str:
        out = []
        while count:
            if self._avail() == 0:              # a string split over a CONTINUE gets a fresh flag byte
                self.i, self.pos = self.i + 1, 0
                wide = bool(self.take(1)[0] & 1)
            width = 2 if wide else 1
            k = min(count, self._avail() // width)
            raw = self.parts[self.i][self.pos:self.pos + k * width]
            self.pos += k * width
            out.append(raw.decode("utf-16-le") if wide else raw.decode("latin-1"))
            count -= k
        return "".join(out)


def _xl_string(p: _Pieces, len_bytes: int = 2) -> str:
    n = struct.unpack("<H" if len_bytes == 2 else "<B", p.take(len_bytes))[0]
    flags = p.take(1)[0]
    runs = struct.unpack("<H", p.take(2))[0] if flags & 8 else 0
    ext = struct.unpack("<I", p.take(4))[0] if flags & 4 else 0
    s = p.chars(n, bool(flags & 1))
    p.take(4 * runs + ext)
    return s


def read_xls(path: Path | str | bytes) -> dict[str, list[list]]:
    data = path if isinstance(path, bytes) else Path(path).read_bytes()
    wb = _cfb_stream(data, ("Workbook", "Book"))
    records: list[tuple[int, int, bytes]] = []
    off = 0
    while off + 4 <= len(wb):
        rid, ln = struct.unpack_from("<HH", wb, off)
        records.append((off, rid, wb[off + 4:off + 4 + ln]))
        off += 4 + ln
    sst: list[str] = []
    sheets: list[tuple[str, int]] = []
    for k, (_, rid, body) in enumerate(records):
        if rid == 0x85:
            pos = struct.unpack_from("<I", body, 0)[0]
            sheets.append((_xl_string(_Pieces([body[6:]]), 1), pos))
        elif rid == 0xFC:
            parts = [body]
            j = k + 1
            while j < len(records) and records[j][1] == 0x3C:
                parts.append(records[j][2])
                j += 1
            p = _Pieces(parts)
            _, unique = struct.unpack("<II", p.take(8))
            sst = [_xl_string(p) for _ in range(unique)]
    by_offset = {o: k for k, (o, _, _) in enumerate(records)}
    out: dict[str, list[list]] = {}
    for name, pos in sheets:
        rows: dict[int, dict[int, object]] = {}
        k = by_offset[pos] + 1
        pending = None
        while k < len(records) and records[k][1] != 0x0A:
            _, rid, b = records[k]
            k += 1
            if rid == 0xFD:
                r, c, _, i = struct.unpack_from("<HHHI", b)
                rows.setdefault(r, {})[c] = sst[i]
            elif rid == 0x203:
                r, c, _, v = struct.unpack_from("<HHHd", b)
                rows.setdefault(r, {})[c] = v
            elif rid == 0x27E:
                r, c, _, v = struct.unpack_from("<HHHI", b)
                rows.setdefault(r, {})[c] = _rk(v)
            elif rid == 0xBD:
                r, c0 = struct.unpack_from("<HH", b)
                n = (len(b) - 6) // 6
                for j in range(n):
                    v = struct.unpack_from("<I", b, 4 + 6 * j + 2)[0]
                    rows.setdefault(r, {})[c0 + j] = _rk(v)
            elif rid == 0x204:
                r, c, _ = struct.unpack_from("<HHH", b)
                rows.setdefault(r, {})[c] = _xl_string(_Pieces([b[6:]]))
            elif rid == 0x205:
                r, c, _, v, is_err = struct.unpack_from("<HHHBB", b)
                rows.setdefault(r, {})[c] = None if is_err else bool(v)
            elif rid == 0x06:
                r, c, _ = struct.unpack_from("<HHH", b)
                res = b[6:14]
                if res[6:8] == b"\xff\xff":
                    pending = (r, c) if res[0] == 0 else None
                    if res[0] == 1:
                        rows.setdefault(r, {})[c] = bool(res[2])
                else:
                    rows.setdefault(r, {})[c] = struct.unpack("<d", res)[0]
            elif rid == 0x207 and pending is not None:
                rows.setdefault(pending[0], {})[pending[1]] = _xl_string(_Pieces([b]))
                pending = None
        out[name] = _dense(rows)
    return out


def read_sheets(path: Path | str) -> dict[str, list[list]]:
    head = Path(path).read_bytes()[:8]
    return read_xls(path) if head == bytes.fromhex("D0CF11E0A1B11AE1") else read_xlsx(path)
