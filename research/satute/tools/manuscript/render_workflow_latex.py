#!/usr/bin/env python3
"""Render Figure 1 from one SVG truth with genuine LaTeX equations.

The R script draws the non-mathematical artwork.  This compositor inserts
path-only LaTeX SVG elements, following the figure workflow used by the
``Is it rotation`` manuscript, and derives the publication PDF and PNG from
that completed SVG.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))

from tools.figures.export import export_svg_pdf, export_svg_png  # noqa: E402
from tools.figures.latex_svg import compile_latex_svg, latex_svg_element  # noqa: E402


WIDTH = 511.0
HEIGHT = 504.0
INK = "#202020"
MID = "#777777"
GREEN = "#58767C"
R_RENDERER = PROJECT / "tools" / "manuscript" / "render_workflow.R"
OUTPUT_DIR = PROJECT / "artifacts" / "releases" / "generated-figures"
OUTPUT_STEM = OUTPUT_DIR / "figure1_satute_weighted_workflow"
MANUSCRIPT_PDF = PROJECT / "manuscript" / "figures" / f"{OUTPUT_STEM.name}.pdf"


@dataclass(frozen=True)
class Equation:
    expression: str
    x: float
    y: float
    font_size: float
    kind: str
    colour: str = INK
    anchor: str = "middle"
    aria_label: str | None = None


# Coordinates reproduce the grid locations in render_workflow.R.  Expressions
# are the manuscript LaTeX, not translations into R plotmath.
EQUATIONS = (
    Equation(r"\partial", 0.255, 0.785, 16.0, "latex-site-pattern"),
    Equation(r"\partial_A", 0.489, 0.704, 11.5, "latex-left-subpattern"),
    Equation(r"\partial_B", 0.882, 0.704, 11.5, "latex-right-subpattern"),
    Equation(r"T_A", 0.600, 0.716, 13.0, "latex-left-subtree"),
    Equation(r"T_B", 0.820, 0.716, 13.0, "latex-right-subtree"),
    Equation(r"L(\partial_A)", 0.642, 0.774, 9.0, "latex-left-likelihood"),
    Equation(r"L(\partial_B)", 0.743, 0.774, 9.0, "latex-right-likelihood"),
    Equation(
        r"C_{\partial i}="
        r"\left(\frac{L(\partial_A)}{P(\partial_A)}\!\cdot h_i\right)"
        r"\left(\frac{L(\partial_B)}{P(\partial_B)}\!\cdot h_i\right)",
        0.696,
        0.505,
        7.7,
        "latex-coherence-coefficient",
    ),
    Equation(
        r"w_i(t)=\exp\{(\lambda_i-\lambda_1)t\}",
        0.660,
        0.322,
        9.0,
        "latex-persistence-weight",
    ),
    Equation(
        r"C_{\partial}^{\lambda}(t)=\sum_{i=1}^{d}w_i(t)C_{\partial i}",
        0.660,
        0.254,
        7.6,
        "latex-weighted-coherence",
    ),
    Equation(
        r"\widehat C_{\lambda}(t)="
        r"\sum_{\partial}\frac{n_{\partial}}{n}C_{\partial}^{\lambda}(t)",
        0.660,
        0.183,
        7.4,
        "latex-alignment-average",
    ),
    Equation(
        r"Z_{\lambda}(t)="
        r"\frac{\sqrt n\,\widehat C_{\lambda}(t)}{\widehat\sigma_{\lambda}(t)}",
        0.660,
        0.121,
        8.4,
        "latex-weighted-statistic",
    ),
    Equation(
        r"H_0:\ \widehat C_{\lambda}(t)"
        r"\mathrel{\dot\sim}\mathcal N\!\left(0,\widehat\sigma_{\lambda}^2(t)/n\right)",
        0.660,
        0.078,
        6.6,
        "latex-null-distribution",
    ),
    Equation(r"z_{\alpha}", 0.3987, 0.126, 7.2, "latex-critical-value", MID),
    Equation(r"\alpha=0.01", 0.4158, 0.178, 7.2, "latex-alpha", GREEN),
    Equation(
        r"Z_{\lambda}(t)>z_{\alpha}",
        0.420,
        0.072,
        10.5,
        "latex-decision-rule",
    ),
)


def require_program(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        raise RuntimeError(f"Figure 1 rendering requires {name}")
    return path


def render_base_svg(svg_path: Path) -> None:
    environment = os.environ.copy()
    environment["SATUTE_LATEX_OVERLAY"] = "1"
    subprocess.run(
        [require_program("Rscript"), str(R_RENDERER)],
        cwd=PROJECT,
        env=environment,
        check=True,
    )
    if not svg_path.is_file():
        raise RuntimeError(f"R renderer did not create {svg_path}")


def equation_element(spec: Equation) -> str:
    centre_x = spec.x * WIDTH
    centre_y = (1.0 - spec.y) * HEIGHT
    compiled = compile_latex_svg(spec.expression, spec.colour)
    # The source document compiles at 10 pt.  Scaling by the nominal font size,
    # rather than by the complete bounding-box height, keeps fractions and
    # display-style sums at the same optical type size as simple expressions.
    height = compiled.height_pt * spec.font_size / 10.0
    return latex_svg_element(
        spec.expression,
        x=centre_x,
        y=centre_y - height / 2.0,
        height=height,
        anchor=spec.anchor,
        colour=spec.colour,
        data_kind=spec.kind,
        aria_label=spec.aria_label,
    )


def insert_equations(svg_path: Path) -> None:
    svg = svg_path.read_text(encoding="utf-8")
    if svg.count("</svg>") != 1:
        raise RuntimeError("Figure 1 base must contain one root SVG")
    overlay = "\n<g id=\"latex-equation-layer\" data-kind=\"latex-equation-layer\">\n"
    overlay += "\n".join(equation_element(spec) for spec in EQUATIONS)
    overlay += "\n</g>\n"
    svg = svg.replace("</svg>", overlay + "</svg>", 1)
    svg_path.write_text(svg, encoding="utf-8")


def export(svg_path: Path, pdf_path: Path, png_path: Path) -> None:
    export_svg_pdf(svg_path, pdf_path)
    export_svg_png(svg_path, png_path, width=2272)


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    svg_path = OUTPUT_STEM.with_suffix(".svg")
    pdf_path = OUTPUT_STEM.with_suffix(".pdf")
    png_path = OUTPUT_STEM.with_suffix(".png")
    render_base_svg(svg_path)
    insert_equations(svg_path)
    export(svg_path, pdf_path, png_path)
    shutil.copy2(pdf_path, MANUSCRIPT_PDF)
    for path in (svg_path, pdf_path, png_path, MANUSCRIPT_PDF):
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
