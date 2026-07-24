#!/usr/bin/env python3
"""Verify that the manuscript contains a complete, readable power-plot atlas."""

from __future__ import annotations

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
    "dominant": "Published dominant coefficient",
    "weighted": "Eigenvalue-weighted coefficient",
    "difference": "Paired formula difference",
}
COMMON_MARKERS = (
    "ML tree",
    "True topology, ML lengths",
    "True tree, fixed lengths",
    "Five-taxon external branch",
    "16-taxon internal branch",
    "taxon_bonferroni",
)


def command_output(*args: str) -> str:
    return subprocess.run(args, check=True, text=True, capture_output=True).stdout


def main() -> int:
    errors: list[str] = []
    for command in ("pdfinfo", "pdftotext"):
        if shutil.which(command) is None:
            errors.append(f"required command is unavailable: {command}")

    supplement_text = SUPPLEMENT.read_text(encoding="utf-8") if SUPPLEMENT.is_file() else ""
    baseline = FIGURE_ROOT / "figure_jc_baseline_power.pdf"
    if not baseline.is_file() or baseline.stat().st_size < 10_000:
        errors.append(f"missing original-paper-style JC benchmark: {baseline}")
    if "figure_jc_baseline_power.pdf" not in supplement_text:
        errors.append("supplement does not reference the original-paper-style JC benchmark")
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
            for data_name in ("head_to_head_summary.tsv", "head_to_head_detail.tsv"):
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
                if path.name not in supplement_text:
                    errors.append(f"supplement does not reference: {path.name}")

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
