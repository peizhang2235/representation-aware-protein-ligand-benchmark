#!/usr/bin/env python3
"""Verify the public archive checksum manifest."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "PUBLIC_MANIFEST_SHA256.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    failures: list[str] = []
    with MANIFEST.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        path = ROOT / row["relative_path"]
        if not path.is_file():
            failures.append(f"missing:{row['relative_path']}")
            continue
        if path.stat().st_size != int(row["bytes"]):
            failures.append(f"size:{row['relative_path']}")
        if sha256(path) != row["sha256"]:
            failures.append(f"sha256:{row['relative_path']}")
    print(f"checked={len(rows)} failures={len(failures)}")
    for item in failures:
        print(item)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
