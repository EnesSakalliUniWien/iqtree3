#!/usr/bin/env python3
"""Render a short LaTeX expression as tightly cropped vector SVG.

Adapted from the ``Is it rotation`` figure toolchain.  The generated SVG uses
vector glyph paths and retains the original expression in
``data-latex-source``.  It is intended for mathematical labels embedded in
figures, not for complete figure layouts.
"""

from __future__ import annotations

import argparse
import hashlib
import html
import re
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path


DEFAULT_INK = "#202020"
PT_TO_MM = 25.4 / 72.27


@dataclass(frozen=True)
class LatexSvg:
    """A compiled, tightly cropped vector expression."""

    svg: str
    width_pt: float
    height_pt: float
    view_box: tuple[float, float, float, float]


def toolchain() -> dict[str, str]:
    resolved = {name: shutil.which(name) for name in ("latex", "dvisvgm")}
    missing = [name for name, path in resolved.items() if path is None]
    if missing:
        raise RuntimeError("LaTeX figure rendering requires " + ", ".join(missing))
    return {name: str(path) for name, path in resolved.items()}


def _dimension_to_pt(value: str) -> float:
    match = re.fullmatch(r"([0-9.+-eE]+)(pt|bp|mm|cm|in)", value.strip())
    if match is None:
        raise ValueError(f"unsupported SVG dimension {value!r}")
    number = float(match.group(1))
    unit = match.group(2)
    return number * {
        "pt": 1.0,
        "bp": 72.27 / 72.0,
        "mm": 1.0 / PT_TO_MM,
        "cm": 10.0 / PT_TO_MM,
        "in": 72.27,
    }[unit]


def _prefix_ids(svg: str, prefix: str) -> str:
    """Prevent dvisvgm glyph identifiers from colliding in a parent SVG."""

    identifiers = re.findall(r'\bid=["\']([^"\']+)["\']', svg)
    for identifier in sorted(set(identifiers), key=len, reverse=True):
        replacement = f"{prefix}-{identifier}"
        svg = re.sub(
            rf'(?<=id=["\']){re.escape(identifier)}(?=["\'])',
            replacement,
            svg,
        )
        svg = svg.replace(f"#{identifier}", f"#{replacement}")
    return svg


@lru_cache(maxsize=128)
def compile_latex_svg(expression: str, colour: str = DEFAULT_INK) -> LatexSvg:
    """Compile one display-style expression to tightly cropped SVG vectors."""

    tools = toolchain()
    hex_colour = colour.removeprefix("#").upper()
    if not re.fullmatch(r"[0-9A-F]{6}", hex_colour):
        raise ValueError("colour must be a six-digit hexadecimal value")
    document = rf"""\documentclass[10pt]{{article}}
\usepackage{{amsmath,amssymb,xcolor}}
\pagestyle{{empty}}
\begin{{document}}
\thispagestyle{{empty}}
\color[HTML]{{{hex_colour}}}$\displaystyle {expression}$
\end{{document}}
"""
    with tempfile.TemporaryDirectory(prefix="satute-figure-latex-") as directory:
        work = Path(directory)
        source = work / "expression.tex"
        source.write_text(document, encoding="utf-8")
        latex = subprocess.run(
            [
                tools["latex"],
                "-interaction=nonstopmode",
                "-halt-on-error",
                "-no-shell-escape",
                source.name,
            ],
            cwd=work,
            capture_output=True,
            text=True,
        )
        if latex.returncode:
            raise RuntimeError(f"LaTeX expression failed:\n{latex.stdout[-3000:]}")
        svg_path = work / "expression.svg"
        convert = subprocess.run(
            [
                tools["dvisvgm"],
                "--exact-bbox",
                "--no-fonts",
                "--precision=6",
                f"--output={svg_path}",
                str(work / "expression.dvi"),
            ],
            cwd=work,
            capture_output=True,
            text=True,
        )
        if convert.returncode:
            raise RuntimeError(f"dvisvgm conversion failed:\n{convert.stderr[-3000:]}")
        raw = svg_path.read_text(encoding="utf-8")

    root = ET.fromstring(raw)
    if "viewBox" not in root.attrib:
        raise ValueError("dvisvgm returned SVG without a viewBox")
    view_box = tuple(float(value) for value in root.attrib["viewBox"].split())
    if len(view_box) != 4:
        raise ValueError("dvisvgm returned an invalid viewBox")
    prefix = "tex-" + hashlib.sha256(f"{expression}|{colour}".encode()).hexdigest()[:12]
    raw = re.sub(r"<\?xml[^>]*>\s*", "", raw, count=1)
    raw = _prefix_ids(raw, prefix)
    return LatexSvg(
        svg=raw,
        width_pt=_dimension_to_pt(root.attrib["width"]),
        height_pt=_dimension_to_pt(root.attrib["height"]),
        view_box=(view_box[0], view_box[1], view_box[2], view_box[3]),
    )


def latex_svg_element(
    expression: str,
    *,
    x: float,
    y: float,
    height: float,
    anchor: str = "middle",
    colour: str = DEFAULT_INK,
    data_kind: str = "latex-equation",
    aria_label: str | None = None,
) -> str:
    """Place a compiled expression by its top edge in parent SVG units."""

    compiled = compile_latex_svg(expression, colour)
    scale = height / compiled.height_pt
    width = compiled.width_pt * scale
    if anchor == "middle":
        left = x - width / 2.0
    elif anchor == "end":
        left = x - width
    elif anchor == "start":
        left = x
    else:
        raise ValueError("anchor must be start, middle or end")
    nested = re.sub(r"\s(?:width|height)=[\"'][^\"']+[\"']", "", compiled.svg)
    semantic_label = aria_label or expression
    element_id = (
        "latex-"
        + re.sub(r"[^a-z0-9-]+", "-", data_kind.lower()).strip("-")
        + "-"
        + hashlib.sha256(expression.encode()).hexdigest()[:8]
    )
    return nested.replace(
        "<svg ",
        (
            f'<svg id="{element_id}" x="{left:.6f}" y="{y:.6f}" '
            f'width="{width:.6f}" height="{height:.6f}" '
            f'data-kind="{html.escape(data_kind, quote=True)}" '
            f'data-latex-source="{html.escape(expression, quote=True)}" '
            f'role="math" aria-label="{html.escape(semantic_label, quote=True)}" '
        ),
        1,
    )


def compile_expression(expression: str, colour: str = DEFAULT_INK) -> str:
    """Compile one display-style expression and return standalone SVG text."""

    compiled = compile_latex_svg(expression, colour)
    svg = compiled.svg.replace(
        "<svg ",
        (
            '<svg data-kind="latex-equation" '
            f'data-latex-source="{html.escape(expression, quote=True)}" '
            f'role="math" aria-label="{html.escape(expression, quote=True)}" '
        ),
        1,
    )
    return '<?xml version="1.0" encoding="UTF-8"?>\n' + svg


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("expression", help="LaTeX math expression without dollar signs")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--colour", default=DEFAULT_INK)
    args = parser.parse_args()
    svg = compile_expression(args.expression, args.colour)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(svg + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
