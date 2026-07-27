#!/usr/bin/env python3
"""Verify that the manuscript contains a complete, readable power-plot atlas."""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
FIGURE_ROOT = PROJECT_ROOT / "manuscript" / "figures"
POWER_ROOT = FIGURE_ROOT / "publication-suite" / "power"
CASE_ROOT = PROJECT_ROOT / "artifacts" / "publication" / "power"
ORIGINAL_ROOT = CASE_ROOT / "original-paper"
SUPPLEMENT = PROJECT_ROOT / "manuscript" / "sections" / "06-supplementary-information.tex"
BASELINE_RESULTS = (
    PROJECT_ROOT / "manuscript" / "sections" / "results" / "22-weighted-simulations.tex"
)
BASELINE_RENDERER = (
    PROJECT_ROOT / "experiments" / "001_paper_reconstruction" / "render_full_baseline.R"
)

CASES = (
    "jc_JC",
    "gtr_GTR_PF06346",
    "gtr_GTR_PF06346_G4",
    "gtr_GTR_PF06346_I_G4",
    "protein_LG_G4",
    "protein_WAG_G4",
    "protein_JTT_G4",
    "protein_Q_PFAM_G4",
    "misspec_JC",
    "misspec_K2P",
    "misspec_F81",
)
PLOT_MARKERS = {
    "dominant": "Published SatuTe statistic",
    "weighted": "Eigenvalue-weighted SatuTe statistic",
    "difference": "Paired formula difference",
}
COMMON_MARKERS = (
    "ML tree",
    "Five-taxon external branch",
    "16-taxon internal branch",
    "taxon-pair Bonferroni",
)
FORBIDDEN_BONFERRONI_MARKERS = (
    "True topology, ML lengths",
    "True tree, fixed lengths",
)
INTEGRATED_FIGURES = {
    "figure_model_misspecification.pdf": (
        CASE_ROOT.parent / "figures" / "model-misspecification",
        (
            "SatuTe under substitution-model misspecification",
            "Bonferroni is used only for the inferred-ML-tree series",
        ),
        ("data/figure_model_misspecification_data.tsv",),
    ),
    "figure_power_curve_summary.pdf": (
        CASE_ROOT.parent / "figures" / "cross-model-power",
        (
            "Cross-model SatuTe power under inferred-tree testing",
            "taxon-pair Bonferroni",
            "curves coincide",
        ),
        (
            "data/figure_power_curve_summary_data.tsv",
            "data/figure_power_curve_summary_claims.tsv",
        ),
    ),
}


def command_output(*args: str) -> str:
    return subprocess.run(args, check=True, text=True, capture_output=True).stdout


def main() -> int:
    errors: list[str] = []
    for command in ("pdfinfo", "pdftotext"):
        if shutil.which(command) is None:
            errors.append(f"required command is unavailable: {command}")

    supplement_text = SUPPLEMENT.read_text(encoding="utf-8") if SUPPLEMENT.is_file() else ""
    baseline_results_text = (
        BASELINE_RESULTS.read_text(encoding="utf-8") if BASELINE_RESULTS.is_file() else ""
    )
    baseline = FIGURE_ROOT / "figure_jc_baseline_power.pdf"
    if not baseline.is_file() or baseline.stat().st_size < 10_000:
        errors.append(f"missing original-paper-style JC benchmark: {baseline}")
    elif not errors:
        baseline_info = command_output("pdfinfo", str(baseline))
        baseline_text = command_output("pdftotext", str(baseline), "-")
        for marker in (
            "Informative instances (%)",
            "Branch length",
            "true tree +",
            "ML-tree",
            "Bonferroni",
        ):
            if marker not in baseline_text:
                errors.append(f"original-style JC benchmark is missing text marker: {marker}")
        size_line = next(
            (line for line in baseline_info.splitlines() if line.startswith("Page size:")),
            "",
        )
        size_match = re.search(r"Page size:\s+([0-9.]+) x ([0-9.]+)", size_line)
        if not size_match or not (
            500 <= float(size_match.group(1)) <= 520
            and 830 <= float(size_match.group(2)) <= 850
        ):
            errors.append(
                "original-style JC benchmark must retain the original portrait page geometry"
            )
    canonical_baseline = CASE_ROOT / "jc_JC" / "figures" / "original-style-baseline.pdf"
    if not canonical_baseline.is_file():
        errors.append(f"missing canonical original-style JC benchmark: {canonical_baseline}")
    elif baseline.is_file() and canonical_baseline.read_bytes() != baseline.read_bytes():
        errors.append("manuscript JC benchmark differs from its canonical case-local figure")
    if "figure_jc_baseline_power.pdf" not in baseline_results_text:
        errors.append("main Results do not reference the original-paper-style JC benchmark")
    renderer_text = (
        BASELINE_RENDERER.read_text(encoding="utf-8") if BASELINE_RENDERER.is_file() else ""
    )
    for marker in (
        "library(ape)",
        "library(ggtree)",
        'layout = "unrooted"',
        "groupOTU(",
    ):
        if marker not in renderer_text:
            errors.append(f"original-style renderer is missing required tree contract: {marker}")
    for reference in ("simulation-main.pdf", "simulation-misspecification.pdf"):
        path = ORIGINAL_ROOT / "figures" / reference
        if not path.is_file():
            errors.append(f"missing indexed original-paper figure: {path}")
    for data_name in (
        "jc_summary.tsv",
        "jc_detail.tsv",
        "misspecification_summary.tsv",
        "misspecification_detail.tsv",
        "archived_methods.tex",
    ):
        path = ORIGINAL_ROOT / "data" / data_name
        if not path.is_file():
            errors.append(f"missing original-comparison data link: {path}")

    if not errors:
        for case in CASES:
            case_dir = CASE_ROOT / case
            for data_name in (
                "head_to_head_summary.tsv",
                "head_to_head_detail.tsv",
                "paired_formula_difference.tsv",
            ):
                path = case_dir / "data" / data_name
                if not path.is_file():
                    errors.append(f"{case} is missing case-local source data: {path}")
            if not (case_dir / "provenance.json").is_file():
                errors.append(f"{case} is missing provenance.json")
            for plot_kind, title_marker in PLOT_MARKERS.items():
                path = POWER_ROOT / f"{case}_{plot_kind}.pdf"
                canonical_path = case_dir / "figures" / f"{plot_kind}.pdf"
                if not canonical_path.is_file() or canonical_path.stat().st_size < 10_000:
                    errors.append(f"missing case-local canonical plot: {canonical_path}")
                if not path.is_file() or path.stat().st_size < 10_000:
                    errors.append(f"missing or undersized plot: {path}")
                    continue
                if canonical_path.is_file() and canonical_path.read_bytes() != path.read_bytes():
                    errors.append(f"manuscript plot differs from canonical plot: {path}")
                try:
                    info = command_output("pdfinfo", str(path))
                    text = command_output("pdftotext", str(path), "-")
                except subprocess.CalledProcessError as exc:
                    errors.append(f"could not inspect {path}: {exc}")
                    continue
                if "Pages:           1" not in info:
                    errors.append(f"plot must be a one-page PDF: {path}")
                for marker in (title_marker, *COMMON_MARKERS):
                    if marker not in text:
                        errors.append(f"{path.name} is missing text marker: {marker}")
                for marker in FORBIDDEN_BONFERRONI_MARKERS:
                    if marker in text:
                        errors.append(
                            f"{path.name} applies taxon-pair Bonferroni to "
                            f"an inapplicable scenario: {marker}"
                        )
                if path.name not in supplement_text:
                    errors.append(f"supplement does not reference: {path.name}")

        results_text = (
            PROJECT_ROOT / "manuscript" / "sections" / "results" / "24-publication-suite.tex"
        ).read_text(encoding="utf-8")
        for filename, (figure_root, markers, data_paths) in INTEGRATED_FIGURES.items():
            canonical_path = figure_root / "figures" / filename
            manuscript_path = FIGURE_ROOT / filename
            if not canonical_path.is_file():
                errors.append(f"missing integrated canonical figure: {canonical_path}")
                continue
            if not manuscript_path.is_file():
                errors.append(f"missing integrated manuscript figure: {manuscript_path}")
                continue
            if canonical_path.read_bytes() != manuscript_path.read_bytes():
                errors.append(f"integrated manuscript figure is stale: {manuscript_path}")
            text = command_output("pdftotext", str(canonical_path), "-")
            for marker in markers:
                if marker not in text:
                    errors.append(f"{filename} is missing text marker: {marker}")
            for data_path in data_paths:
                path = figure_root / data_path
                if not path.is_file() or path.stat().st_size == 0:
                    errors.append(f"{filename} is missing recalculation data: {path}")
            if filename not in results_text:
                errors.append(f"main Results do not reference integrated figure: {filename}")

    if errors:
        print("Power-plot atlas audit failed:", file=sys.stderr)
        for error in errors:
            print(f"  - {error}", file=sys.stderr)
        return 1

    print(
        "Power-plot atlas audit passed: "
        f"{len(CASES)} cases × {len(PLOT_MARKERS)} views, plus the original-style JC benchmark."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
