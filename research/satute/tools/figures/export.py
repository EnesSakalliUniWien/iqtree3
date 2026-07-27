#!/usr/bin/env python3
"""Publication-safe SVG export with embedded ordinary-text fonts."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from pathlib import Path


def _prepare_cairo_runtime() -> None:
    homebrew_library = Path("/opt/homebrew/lib")
    if not homebrew_library.is_dir():
        return
    paths = [
        item
        for item in os.environ.get("DYLD_FALLBACK_LIBRARY_PATH", "").split(":")
        if item
    ]
    if str(homebrew_library) not in paths:
        os.environ["DYLD_FALLBACK_LIBRARY_PATH"] = ":".join(
            [str(homebrew_library), *paths]
        )


def export_svg_pdf(svg_path: Path, pdf_path: Path) -> None:
    """Export SVG to vector PDF while embedding ordinary text fonts."""

    _prepare_cairo_runtime()
    try:
        import cairosvg
    except (ImportError, OSError) as error:
        raise RuntimeError("SVG PDF export requires CairoSVG and Cairo") from error
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    cairosvg.svg2pdf(url=str(svg_path), write_to=str(pdf_path))
    qpdf = shutil.which("qpdf")
    if qpdf is None:
        raise RuntimeError("SVG PDF export requires qpdf")
    subprocess.run(
        [qpdf, "--check", str(pdf_path)],
        check=True,
        capture_output=True,
        text=True,
    )


def export_svg_png(
    svg_path: Path,
    png_path: Path,
    *,
    scale: float | None = None,
    width: int | None = None,
) -> None:
    """Rasterize the same SVG truth to PNG."""

    converter = shutil.which("rsvg-convert")
    if converter is None:
        raise RuntimeError("SVG PNG export requires rsvg-convert")
    command = [converter, "--format=png"]
    if scale is not None:
        command.extend(["--zoom", str(scale)])
    if width is not None:
        command.extend(["--width", str(width)])
    png_path.parent.mkdir(parents=True, exist_ok=True)
    command.extend([f"--output={png_path}", str(svg_path)])
    subprocess.run(command, check=True, capture_output=True, text=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("svg", type=Path)
    parser.add_argument("--pdf", type=Path, required=True)
    parser.add_argument("--png", type=Path)
    parser.add_argument("--png-scale", type=float)
    parser.add_argument("--png-width", type=int)
    args = parser.parse_args()
    export_svg_pdf(args.svg, args.pdf)
    if args.png is not None:
        export_svg_png(
            args.svg,
            args.png,
            scale=args.png_scale,
            width=args.png_width,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
