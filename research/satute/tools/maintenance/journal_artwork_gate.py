#!/usr/bin/env python3
"""Audit registered PDF figures for basic journal-production properties.

Adapted from the journal-artwork gate in the ``Is it rotation`` manuscript.
The SatuTe figures are PDF outputs from R and Python, so this gate checks the
portable PDF contract rather than SVG-specific illustration metadata.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[2]
REGISTRY = PROJECT / "manuscript" / "figures" / "registry.json"
MAXIMUM_BYTES = 50 * 1024 * 1024
MINIMUM_DIMENSION_PT = 100.0
MAXIMUM_DIMENSION_PT = 2000.0
FONT_REQUIRED_ASSETS = {
    "manuscript/figures/figure1_satute_weighted_workflow.pdf",
}


def output(*command: str) -> str:
    return subprocess.run(
        command,
        cwd=PROJECT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout


def font_embedding_failures(path: Path, *, require_fonts: bool) -> list[str]:
    lines = output("pdffonts", str(path)).splitlines()[2:]
    failures: list[str] = []
    if require_fonts and not lines:
        failures.append(f"{path}: ordinary figure text was not retained as embedded fonts")
    for line in lines:
        fields = line.split()
        if len(fields) >= 6 and fields[-5].lower() != "yes":
            failures.append(f"{path}: font is not embedded: {fields[0]}")
        if "Type 3" in line:
            failures.append(f"{path}: Type 3 font is not publication-safe: {fields[0]}")
    return failures


def main() -> int:
    failures: list[str] = []
    for command in ("pdfinfo", "pdffonts", "pdftotext"):
        if shutil.which(command) is None:
            failures.append(f"required production tool is unavailable: {command}")
    if failures:
        for failure in failures:
            print(f"  - {failure}", file=sys.stderr)
        return 1

    registry = json.loads(REGISTRY.read_text(encoding="utf-8"))
    assets = sorted(
        {
            str(asset)
            for entry in registry["figures"]
            for asset in entry["assets"]
        }
    )
    for relative in assets:
        path = PROJECT / relative
        if not path.is_file():
            failures.append(f"missing registered PDF: {relative}")
            continue
        if path.suffix.lower() != ".pdf":
            failures.append(f"registered manuscript asset is not PDF: {relative}")
            continue
        if path.stat().st_size == 0 or path.stat().st_size > MAXIMUM_BYTES:
            failures.append(f"invalid file size for {relative}: {path.stat().st_size} bytes")
        try:
            info = output("pdfinfo", str(path))
            text = output("pdftotext", str(path), "-")
            failures.extend(
                font_embedding_failures(
                    path,
                    require_fonts=relative in FONT_REQUIRED_ASSETS,
                )
            )
        except subprocess.CalledProcessError as error:
            failures.append(f"could not inspect {relative}: {error}")
            continue

        pages = re.search(r"^Pages:\s+(\d+)", info, re.MULTILINE)
        if not pages or int(pages.group(1)) != 1:
            failures.append(f"{relative}: figure PDF must contain exactly one page")
        size = re.search(
            r"^Page size:\s+([0-9.]+) x ([0-9.]+) pts",
            info,
            re.MULTILINE,
        )
        if not size:
            failures.append(f"{relative}: page dimensions are unavailable")
        else:
            width, height = map(float, size.groups())
            if not (
                MINIMUM_DIMENSION_PT <= width <= MAXIMUM_DIMENSION_PT
                and MINIMUM_DIMENSION_PT <= height <= MAXIMUM_DIMENSION_PT
            ):
                failures.append(f"{relative}: implausible page dimensions {width:g} x {height:g} pt")
        if len(re.findall(r"[A-Za-z0-9]+", text)) < 2:
            failures.append(f"{relative}: PDF has insufficient extractable text")

    if failures:
        print("Journal artwork audit failed:", file=sys.stderr)
        for failure in failures:
            print(f"  - {failure}", file=sys.stderr)
        return 1
    print(f"Journal artwork audit passed: {len(assets)} registered PDF assets.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
