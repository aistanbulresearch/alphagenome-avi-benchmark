"""Download public data files and record where they came from.

Each download is appended to data/downloads.json with the URL, the HTTP Last-Modified header, the byte
count, the SHA-256 of the saved file and the UTC start and end times. An existing file is never
overwritten.

Usage:
    python src/download_public.py <url> [<url> ...]
    python src/download_public.py --name <relative/path/in/raw/folder> <url>
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

from paths import RAW as RAW_DIR, ROOT as STUDY_ROOT

LEDGER = STUDY_ROOT / "data" / "downloads.json"
CHUNK = 1 << 20


def _now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def target_path(url: str, name: str | None) -> Path:
    """<raw folder>/<name>, or <raw folder>/<last URL segment>; a name must stay inside the raw folder."""
    if name is None:
        return RAW_DIR / url.rstrip("/").rsplit("/", 1)[-1]
    rel = Path(name)
    target = RAW_DIR / rel
    if (rel.is_absolute() or rel.drive or rel.root or ".." in rel.parts or not rel.parts
            or RAW_DIR.resolve() not in target.resolve().parents):
        raise ValueError(f"name must be a relative path inside the raw folder: {name!r}")
    return target


def download(url: str, name: str | None = None) -> dict:
    """Download one URL into the raw folder and return its record."""
    target = target_path(url, name)
    if target.exists():
        raise FileExistsError(f"{target} already exists; refusing to overwrite")
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(target.suffix + ".part")
    started = _now()
    sha = hashlib.sha256()
    size = 0
    request = urllib.request.Request(url, headers={"User-Agent": "alphagenome-avi-benchmark/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response, open(partial, "wb") as out:
        status = response.status
        last_modified = response.headers.get("Last-Modified")
        declared = response.headers.get("Content-Length")
        while True:
            block = response.read(CHUNK)
            if not block:
                break
            out.write(block)
            sha.update(block)
            size += len(block)
    if declared is not None and int(declared) != size:
        raise IOError(f"{url}: received {size} bytes, server declared {declared}")
    partial.rename(target)
    return {
        "url": url,
        "file": "data/raw/" + str(target.relative_to(RAW_DIR)).replace("\\", "/"),
        "http_status": status,
        "last_modified": last_modified,
        "bytes": size,
        "sha256": sha.hexdigest(),
        "started_utc": started,
        "finished_utc": _now(),
    }


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--name", help="target path inside the raw folder (one URL only)")
    ap.add_argument("urls", nargs="+")
    args = ap.parse_args(argv)
    if args.name is not None and len(args.urls) != 1:
        raise ValueError("--name takes exactly one URL")
    ledger = json.loads(LEDGER.read_text(encoding="utf-8")) if LEDGER.exists() else []
    for url in args.urls:
        record = download(url, args.name)
        ledger.append(record)
        LEDGER.parent.mkdir(parents=True, exist_ok=True)
        LEDGER.write_text(json.dumps(ledger, indent=2) + "\n", encoding="utf-8", newline="\n")
        print(json.dumps(record))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
