#!/usr/bin/env python3
"""Rebuild compact visual-QC sheets from the current v1.7 PNG figures."""

from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path

from PIL import Image


PACKAGE = Path(__file__).resolve().parents[1]
QC = PACKAGE / "qc"


def figure_number(path: Path) -> int:
    match = re.search(r"Figure_(\d+)", path.name)
    if match is None:
        raise ValueError(f"Cannot parse figure number: {path.name}")
    return int(match.group(1))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_sheet(
    source_directory: Path,
    pattern: str,
    expected_count: int,
    columns: int,
    thumbnail_width: int,
    gap: int,
    destination: Path,
) -> dict[str, object]:
    sources = sorted(source_directory.glob(pattern), key=figure_number)
    if len(sources) != expected_count:
        raise RuntimeError(
            f"Expected {expected_count} PNGs matching {pattern}, found {len(sources)}"
        )

    thumbnails: list[Image.Image] = []
    for source in sources:
        with Image.open(source) as image:
            height = round(image.height * thumbnail_width / image.width)
            thumbnail = image.convert("RGB").resize(
                (thumbnail_width, height), Image.Resampling.LANCZOS
            )
        thumbnails.append(thumbnail)

    rows = math.ceil(len(thumbnails) / columns)
    row_heights = [
        max(
            thumbnail.height
            for thumbnail in thumbnails[row * columns : (row + 1) * columns]
        )
        for row in range(rows)
    ]
    canvas_width = columns * thumbnail_width + (columns + 1) * gap
    canvas_height = sum(row_heights) + (rows + 1) * gap
    canvas = Image.new("RGB", (canvas_width, canvas_height), "white")

    y = gap
    for row, row_height in enumerate(row_heights):
        for column, thumbnail in enumerate(
            thumbnails[row * columns : (row + 1) * columns]
        ):
            x = gap + column * (thumbnail_width + gap)
            centered_y = y + (row_height - thumbnail.height) // 2
            canvas.paste(thumbnail, (x, centered_y))
        y += row_height + gap

    destination.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(destination, format="PNG", optimize=True, dpi=(150, 150))
    return {
        "output": str(destination.relative_to(PACKAGE)),
        "figure_count": len(sources),
        "width": canvas.width,
        "height": canvas.height,
        "sha256": sha256(destination),
    }


def main() -> None:
    QC.mkdir(parents=True, exist_ok=True)
    result = {
        "main": build_sheet(
            PACKAGE / "figures" / "main",
            "Figure_*_v1_7_*.png",
            expected_count=6,
            columns=2,
            thumbnail_width=1200,
            gap=16,
            destination=QC / "V1_7_MAIN_FIGURE_CONTACT_SHEET.png",
        ),
        "supplementary": build_sheet(
            PACKAGE / "figures" / "supplementary",
            "Supplementary_Figure_*_v1_7_*.png",
            expected_count=11,
            columns=3,
            thumbnail_width=900,
            gap=15,
            destination=QC / "V1_7_SUPPLEMENTARY_FIGURE_CONTACT_SHEET.png",
        ),
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
