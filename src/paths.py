"""Project paths and file hashing.

Raw downloads are not stored in the repository. They are read from `data/raw` unless the environment
variable AVI_RAW_DIR names another folder.
"""
from __future__ import annotations

import hashlib
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW = Path(os.environ.get("AVI_RAW_DIR") or ROOT / "data" / "raw")
DERIVED = ROOT / "data" / "derived"
EXPOSURE = DERIVED / "exposure"


def sha256_file(path: Path) -> str:
    sha = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            sha.update(chunk)
    return sha.hexdigest()


__all__ = ["ROOT", "RAW", "DERIVED", "EXPOSURE", "sha256_file"]
