"""Spreadsheet readers (standard library only): .xlsx and .xls (BIFF8 inside an OLE2 compound file).
Each returns {sheet name: list of rows}, a row being a list of cell values (str, float, bool, None)."""
from __future__ import annotations

import io
import struct
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------------------------------
_NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
       "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
       "p": "http://schemas.openxmlformats.org/package/2006/relationships"}


def _col_index(ref: str) -> int:
    n = 0
    for ch in ref:
        if not ch.isalpha():
            break
        n = n * 26 + (ord(ch.upper()) - 64)
    return n - 1


def _dense(rows: dict[int, dict[int, object]]) -> list[list]:
    rows = {r: c for r, c in rows.items() if c}
    if not rows:
        return []
    n_cols = max(max(c) for c in rows.values()) + 1
    grid = [[None] * n_cols for _ in range(max(rows) + 1)]
    for r, cells in rows.items():
        for c, v in cells.items():
            grid[r][c] = v
    return grid


def read_xlsx(source: Path | str | bytes) -> dict[str, list[list]]:
    data = source if isinstance(source, bytes) else Path(source).read_bytes()
    m = f"{{{_NS['m']}}}"
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        shared: list[str] = []
        if "xl/sharedStrings.xml" in zf.namelist():
            for si in ET.fromstring(zf.read("xl/sharedStrings.xml")).findall("m:si", _NS):
                shared.append("".join(t.text or "" for t in si.iter(m + "t")))
        rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        target = {r.get("Id"): r.get("Target").lstrip("/") for r in rels.findall("p:Relationship", _NS)}
        out: dict[str, list[list]] = {}
        for sheet in ET.fromstring(zf.read("xl/workbook.xml")).find("m:sheets", _NS).findall("m:sheet", _NS):
            t = target[sheet.get(f"{{{_NS['r']}}}id")]
            rows: dict[int, dict[int, object]] = {}
            for row in ET.fromstring(zf.read(t if t.startswith("xl/") else "xl/" + t)).iter(m + "row"):
                cells = rows.setdefault(int(row.get("r")) - 1, {})
                for c in row.findall("m:c", _NS):
                    typ, v = c.get("t"), c.find("m:v", _NS)
                    if typ == "inlineStr":
                        val = "".join(x.text or "" for x in c.iter(m + "t"))
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
                    cells[_col_index(c.get("r"))] = val
            out[sheet.get("name")] = _dense(rows)
        return out


def _cfb_stream(data: bytes, wanted: tuple[str, ...]) -> bytes:
    """Return the first stream named in `wanted` from an OLE2 compound file."""
    if data[:8] != bytes.fromhex("D0CF11E0A1B11AE1"):
        raise ValueError("not an OLE2 compound file")
    sec = 1 << struct.unpack_from("<H", data, 30)[0]
    mini = 1 << struct.unpack_from("<H", data, 32)[0]
    n_fat, dir_start = struct.unpack_from("<II", data, 44)
    cutoff, minifat_start, _, difat_start, n_difat = struct.unpack_from("<IIIII", data, 56)

    def sector(i: int) -> bytes:
        return data[(i + 1) * sec:(i + 2) * sec]

    def ints(b: bytes) -> tuple:
        return struct.unpack(f"<{len(b) // 4}I", b)

    difat = list(struct.unpack_from("<109I", data, 76))
    nxt = difat_start
    for _ in range(n_difat):
        block = ints(sector(nxt))
        difat.extend(block[:-1])
        nxt = block[-1]
    fat = [x for s in difat[:n_fat] for x in ints(sector(s))]

    def chain(start: int, table) -> list[int]:
        out, s = [], start
        while s < 0xFFFFFFFA:
            out.append(s)
            s = table[s]
        return out

    directory = b"".join(sector(s) for s in chain(dir_start, fat))
    entries = []
    for off in range(0, len(directory), 128):
        e = directory[off:off + 128]
        n = struct.unpack_from("<H", e, 64)[0]
        start, size = struct.unpack_from("<II", e, 116)
        entries.append((e[:max(n - 2, 0)].decode("utf-16-le", "replace"), e[66], start, size))
    for name, typ, start, size in entries:
        if typ != 2 or name not in wanted:
            continue
        if size >= cutoff:
            return b"".join(sector(s) for s in chain(start, fat))[:size]
        ministream = b"".join(sector(s) for s in chain(entries[0][2], fat))[:entries[0][3]]
        minifat = [x for s in chain(minifat_start, fat) for x in ints(sector(s))]
        return b"".join(ministream[s * mini:(s + 1) * mini] for s in chain(start, minifat))[:size]
    raise ValueError(f"no stream named {wanted} in the compound file")


def _rk(value: int) -> float:
    if value & 2:
        num = float((value >> 2) - (1 << 30) if value & 0x80000000 else value >> 2)
    else:
        num = struct.unpack("<d", struct.pack("<Q", (value & 0xFFFFFFFC) << 32))[0]
    return num / 100 if value & 1 else num


class _Pieces:
    """Byte reader over a record body and its CONTINUE records (shared string table)."""

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
            self.pos, n = self.pos + k, n - k
        return out

    def chars(self, count: int, wide: bool) -> str:
        out = []
        while count:
            if self._avail() == 0:              # a string split over CONTINUE gets a fresh flag byte
                self.i, self.pos = self.i + 1, 0
                wide = bool(self.take(1)[0] & 1)
            width = 2 if wide else 1
            k = min(count, self._avail() // width)
            raw = self.parts[self.i][self.pos:self.pos + k * width]
            self.pos, count = self.pos + k * width, count - k
            out.append(raw.decode("utf-16-le") if wide else raw.decode("latin-1"))
        return "".join(out)


def _xl_string(p: _Pieces, len_bytes: int = 2) -> str:
    n = struct.unpack("<H" if len_bytes == 2 else "<B", p.take(len_bytes))[0]
    flags = p.take(1)[0]
    runs = struct.unpack("<H", p.take(2))[0] if flags & 8 else 0
    ext = struct.unpack("<I", p.take(4))[0] if flags & 4 else 0
    s = p.chars(n, bool(flags & 1))
    p.take(4 * runs + ext)
    return s


def parse_biff(wb: bytes) -> dict[str, list[list]]:
    """Cells of every sheet of a BIFF8 workbook stream."""
    records, off = [], 0
    while off + 4 <= len(wb):
        rid, ln = struct.unpack_from("<HH", wb, off)
        records.append((off, rid, wb[off + 4:off + 4 + ln]))
        off += 4 + ln
    sst: list[str] = []
    sheets: list[tuple[str, int]] = []
    for k, (_, rid, body) in enumerate(records):
        if rid == 0x85:
            sheets.append((_xl_string(_Pieces([body[6:]]), 1), struct.unpack_from("<I", body)[0]))
        elif rid == 0xFC:
            parts, j = [body], k + 1
            while j < len(records) and records[j][1] == 0x3C:
                parts.append(records[j][2])
                j += 1
            p = _Pieces(parts)
            unique = struct.unpack("<II", p.take(8))[1]
            sst = [_xl_string(p) for _ in range(unique)]
    by_offset = {o: k for k, (o, _, _) in enumerate(records)}
    out: dict[str, list[list]] = {}
    for name, pos in sheets:
        rows: dict[int, dict[int, object]] = {}
        k, pending = by_offset[pos] + 1, None
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
                for j in range((len(b) - 6) // 6):
                    rows.setdefault(r, {})[c0 + j] = _rk(struct.unpack_from("<I", b, 6 + 6 * j)[0])
            elif rid == 0x204:
                r, c, _ = struct.unpack_from("<HHH", b)
                rows.setdefault(r, {})[c] = _xl_string(_Pieces([b[6:]]))
            elif rid == 0x205:
                r, c, _, v, is_err = struct.unpack_from("<HHHBB", b)
                rows.setdefault(r, {})[c] = None if is_err else bool(v)
            elif rid == 0x06:
                r, c, _ = struct.unpack_from("<HHH", b)
                res = b[6:14]
                pending = None
                if res[6:8] != b"\xff\xff":
                    rows.setdefault(r, {})[c] = struct.unpack("<d", res)[0]
                elif res[0] == 0:
                    pending = (r, c)
                elif res[0] == 1:
                    rows.setdefault(r, {})[c] = bool(res[2])
            elif rid == 0x207 and pending is not None:
                rows.setdefault(pending[0], {})[pending[1]] = _xl_string(_Pieces([b]))
                pending = None
        out[name] = _dense(rows)
    return out


def read_xls(source: Path | str | bytes) -> dict[str, list[list]]:
    data = source if isinstance(source, bytes) else Path(source).read_bytes()
    return parse_biff(_cfb_stream(data, ("Workbook", "Book")))


def read_sheets(path: Path | str) -> dict[str, list[list]]:
    data = Path(path).read_bytes()
    return read_xls(data) if data[:8] == bytes.fromhex("D0CF11E0A1B11AE1") else read_xlsx(data)


def sheet_frame(rows: list[list], header_row: int, keep_only=None) -> pd.DataFrame:
    """Rows below `header_row` as a frame. With `keep_only` (column names), every other column is
    discarded on reading, so prediction columns of a source table are never loaded."""
    header = rows[header_row]
    wanted = None if keep_only is None else {str(k).strip() for k in keep_only}
    keep = [i for i, h in enumerate(header)
            if h is not None and (wanted is None or str(h).strip() in wanted)]
    body = [[r[i] if i < len(r) else None for i in keep] for r in rows[header_row + 1:]]
    body = [r for r in body if any(v is not None for v in r)]
    return pd.DataFrame(body, columns=[str(header[i]).strip() for i in keep])
