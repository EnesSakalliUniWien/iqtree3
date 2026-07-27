#!/usr/bin/env python3
"""Verify manuscript equations and their binding to explanatory artwork."""

from __future__ import annotations

import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[2]
MANUSCRIPT = PROJECT / "manuscript"
WORKFLOW_SOURCE = PROJECT / "tools" / "manuscript" / "render_workflow_latex.py"
WORKFLOW_SVG = (
    PROJECT
    / "artifacts"
    / "releases"
    / "generated-figures"
    / "figure1_satute_weighted_workflow.svg"
)
WORKFLOW_PDF = MANUSCRIPT / "figures" / "figure1_satute_weighted_workflow.pdf"
APPROVED_FIGURE_FONTS = {"Helvetica", "Courier"}
NUMBERED_DISPLAY_PATTERN = re.compile(
    r"\\begin\{(equation|align|gather|multline)\}(.*?)\\end\{\1\}",
    re.DOTALL,
)
LABEL_PATTERN = re.compile(r"\\label\{(eq-[^}]+)\}")
REFERENCE_PATTERN = re.compile(r"\\(?:[cC]ref|eqref|ref)\{([^}]+)\}")
INPUT_PATTERN = re.compile(r"\\input\{([^}]+)\}")


def active_sources(root: Path) -> list[tuple[Path, str]]:
    output: list[tuple[Path, str]] = []
    seen: set[Path] = set()

    def visit(path: Path) -> None:
        resolved = path.resolve()
        if resolved in seen:
            return
        seen.add(resolved)
        text = resolved.read_text(encoding="utf-8")
        output.append((resolved, text))
        for raw in INPUT_PATTERN.findall(text):
            child = MANUSCRIPT / raw
            if child.suffix != ".tex":
                child = child.with_suffix(".tex")
            visit(child)

    visit(root)
    return output


def normalized(text: str) -> str:
    return re.sub(r"\s+", "", text)


def main() -> int:
    failures: list[str] = []
    sys.path.insert(0, str(PROJECT))
    from tools.manuscript.render_workflow_latex import EQUATIONS

    sources = active_sources(MANUSCRIPT / "main.tex")
    sources.extend(active_sources(MANUSCRIPT / "supplementary.tex"))
    all_text = "\n".join(text for _, text in sources)

    labels: dict[str, Path] = {}
    equation_count = 0
    for path, text in sources:
        for environment, body in NUMBERED_DISPLAY_PATTERN.findall(text):
            equation_count += 1
            found = LABEL_PATTERN.findall(body)
            if len(found) != 1:
                failures.append(
                    f"{path.relative_to(PROJECT)}: numbered {environment} display "
                    "must contain exactly one eq-* label"
                )
                continue
            label = found[0]
            if label in labels:
                failures.append(
                    f"duplicate equation label {label}: "
                    f"{labels[label].relative_to(PROJECT)} and {path.relative_to(PROJECT)}"
                )
            labels[label] = path

    references: set[str] = set()
    for group in REFERENCE_PATTERN.findall(all_text):
        references.update(item.strip() for item in group.split(","))
    unknown = sorted(reference for reference in references if reference.startswith("eq-") and reference not in labels)
    if unknown:
        failures.append(f"references to unknown equations: {unknown}")
    unreferenced = sorted(set(labels) - references)
    if unreferenced:
        failures.append(f"numbered equations are never referenced: {unreferenced}")

    aux_path = MANUSCRIPT / "build" / "main.aux"
    log_path = MANUSCRIPT / "build" / "main.log"
    if not aux_path.is_file() or not log_path.is_file():
        failures.append("compiled manuscript metadata is missing; run make manuscript")
    else:
        aux = aux_path.read_text(encoding="utf-8", errors="replace")
        log = log_path.read_text(encoding="utf-8", errors="replace")
        for label in labels:
            if f"{{{label}}}" not in aux:
                failures.append(f"compiled manuscript does not contain equation label {label}")
        for warning in ("undefined references", "multiply defined"):
            if warning in log.lower():
                failures.append(f"compiled manuscript reports {warning}")

    construction = next(
        (
            text
            for path, text in sources
            if path.name == "20-statistical-construction.tex"
        ),
        "",
    )
    formula_bindings = (
        (
            r"w_i(t)=\exp\{(\lambda_i-\lambda_1)t\}",
            r"w_i(t)=\exp\{(\lambda_i-\lambda_1)t\}",
        ),
        (
            r"C_{\partial}^{\lambda}(t)",
            r"C_{\partial}^{\lambda}(t)=\sum_{i=1}^{d}w_i(t)C_{\partial i}",
        ),
        (
            r"Z_{\lambda}(t)",
            r"Z_{\lambda}(t)=",
        ),
    )
    workflow_equations = {normalized(spec.expression) for spec in EQUATIONS}
    for manuscript_fragment, figure_fragment in formula_bindings:
        if normalized(manuscript_fragment) not in normalized(construction):
            failures.append(f"manuscript equation binding is missing: {manuscript_fragment}")
        if not any(normalized(figure_fragment) in expression for expression in workflow_equations):
            failures.append(f"Figure 1 equation binding is missing: {figure_fragment}")

    if not WORKFLOW_SVG.is_file() or not WORKFLOW_PDF.is_file():
        failures.append("final Figure 1 SVG/PDF is missing; run make workflow-figure")
    else:
        workflow_root = ET.parse(WORKFLOW_SVG).getroot()
        text_nodes = [
            element
            for element in workflow_root.iter()
            if element.tag.rsplit("}", 1)[-1] == "text"
        ]
        if not text_nodes:
            failures.append("Figure 1 ordinary labels are not retained as SVG text")
        for element in text_nodes:
            style = element.attrib.get("style", "")
            family = re.search(r'font-family:\s*["\']?([^;"\']+)', style)
            if family is None:
                failures.append(
                    f"Figure 1 text lacks an explicit font family: {element.text!r}"
                )
            elif family.group(1).strip() not in APPROVED_FIGURE_FONTS:
                failures.append(
                    f"Figure 1 uses an unapproved font "
                    f"{family.group(1).strip()!r}: {element.text!r}"
                )
        rendered = {
            element.attrib["data-latex-source"]: element
            for element in workflow_root.iter()
            if "data-latex-source" in element.attrib
        }
        expected = {spec.expression for spec in EQUATIONS}
        missing = sorted(expected - set(rendered))
        unexpected = sorted(set(rendered) - expected)
        if missing:
            failures.append(f"Figure 1 is missing LaTeX expressions: {missing}")
        if unexpected:
            failures.append(f"Figure 1 contains unregistered LaTeX expressions: {unexpected}")
        if len(rendered) != len(EQUATIONS):
            failures.append(
                f"Figure 1 contains {len(rendered)} LaTeX elements; "
                f"expected {len(EQUATIONS)}"
            )
        for expression, element in rendered.items():
            descendants = list(element.iter())
            if not any(node.tag.rsplit("}", 1)[-1] == "path" for node in descendants):
                failures.append(f"Figure 1 equation lacks vector paths: {expression}")
            if any(node.tag.rsplit("}", 1)[-1] == "text" for node in descendants):
                failures.append(f"Figure 1 equation contains a text node: {expression}")

    try:
        from tools.figures.latex_svg import compile_expression

        smoke_expression = r"Z_{\lambda}(t)=\frac{\sqrt{n}\hat C_{\lambda}(t)}{\hat\sigma_{\lambda}(t)}"
        svg = compile_expression(smoke_expression)
        root = ET.fromstring(svg)
        if root.attrib.get("data-latex-source") != smoke_expression:
            failures.append("LaTeX SVG smoke output did not retain its source expression")
        if not any(element.tag.rsplit("}", 1)[-1] == "path" for element in root.iter()):
            failures.append("LaTeX SVG smoke output does not contain vector glyph paths")
        if any(element.tag.rsplit("}", 1)[-1] == "text" for element in root.iter()):
            failures.append("LaTeX SVG smoke output contains renderer-dependent text nodes")
    except Exception as error:  # toolchain diagnostics need the original error
        failures.append(f"LaTeX-to-SVG smoke test failed: {error}")

    if failures:
        print("Equation integration audit failed:", file=sys.stderr)
        for failure in failures:
            print(f"  - {failure}", file=sys.stderr)
        return 1
    print(
        "Equation integration audit passed: "
        f"{equation_count} numbered equations, {len(labels)} unique labels, "
        f"{len(references & set(labels))} referenced labels, 3 manuscript bindings, "
        f"{len(EQUATIONS)} final Figure 1 LaTeX elements, and approved embedded-font sources."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
