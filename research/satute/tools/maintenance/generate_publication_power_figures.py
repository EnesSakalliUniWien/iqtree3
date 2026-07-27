#!/usr/bin/env python3
"""Build and publish every manuscript power plot from one canonical case map."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from dataclasses import asdict, dataclass
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[2]
RENDERER = PROJECT / "experiments" / "002_relative_weighting" / "render_formula_comparison.R"
BASELINE_RENDERER = (
    PROJECT / "experiments" / "001_paper_reconstruction" / "render_full_baseline.R"
)
MISSPECIFICATION_RENDERER = (
    PROJECT
    / "experiments"
    / "002_relative_weighting"
    / "render_misspecification_summary.py"
)
CROSS_MODEL_RENDERER = (
    PROJECT / "experiments" / "002_relative_weighting" / "render_power_summary.py"
)
INTEGRATED_ROOT = PROJECT / "artifacts" / "publication" / "figures"
CASE_ROOT = PROJECT / "artifacts" / "publication" / "power"
PUBLICATION_ROOT = PROJECT / "manuscript" / "figures" / "publication-suite" / "power"


@dataclass(frozen=True)
class Case:
    id: str
    experiment: str
    simulation_model: str
    evaluation_model: str
    dataset: str


JC_DATASET = "artifacts/cluster/lisc/paper-fig2-jc-reps1000-v4-217bd244e77-20260715/merged"
GTR_DATASET = "artifacts/cluster/lisc/paper-matched-gtr-rate-invar-reps1000-ff7c9034b00-20260715-postprocessed-7b264fe9c5ef-v2/merged"
PROTEIN_DATASET = "artifacts/cluster/lisc/paper-matched-protein-g4-reps1000-ff7c9034b00-20260715-postprocessed-7b264fe9c5ef/merged"
MISSPEC_DATASET = "artifacts/cluster/lisc/paper-gtr-misspec-reps1000-ff7c9034b00-20260715-postprocessed-7b264fe9c5ef-v5/merged"

# This is the sole publication power-figure case map. Add or remove cases here;
# the case directories, data links, rendered figures, and manuscript copies all
# follow from this declaration.
CASES = (
    Case("jc_JC", "001_paper_reconstruction", "JC", "JC", JC_DATASET),
    Case("gtr_GTR_PF06346", "002_relative_weighting", "GTR_PF06346", "GTR_PF06346", GTR_DATASET),
    Case("gtr_GTR_PF06346_G4", "002_relative_weighting", "GTR_PF06346_G4", "GTR_PF06346_G4", GTR_DATASET),
    Case("gtr_GTR_PF06346_I_G4", "002_relative_weighting", "GTR_PF06346_I_G4", "GTR_PF06346_I_G4", GTR_DATASET),
    Case("protein_LG_G4", "002_relative_weighting", "LG_G4", "LG_G4", PROTEIN_DATASET),
    Case("protein_WAG_G4", "002_relative_weighting", "WAG_G4", "WAG_G4", PROTEIN_DATASET),
    Case("protein_JTT_G4", "002_relative_weighting", "JTT_G4", "JTT_G4", PROTEIN_DATASET),
    Case("protein_Q_PFAM_G4", "002_relative_weighting", "Q.PFAM_G4", "Q.PFAM_G4", PROTEIN_DATASET),
    Case("misspec_JC", "002_relative_weighting", "GTR_PF06346", "JC", MISSPEC_DATASET),
    Case("misspec_K2P", "002_relative_weighting", "GTR_PF06346", "K2P", MISSPEC_DATASET),
    Case("misspec_F81", "002_relative_weighting", "GTR_PF06346", "F81", MISSPEC_DATASET),
)
PLOT_FILES = {
    "dominant": "figure_dominant_2d_*.pdf",
    "weighted": "figure_eigenvalue_weighted_2d_*.pdf",
    "difference": "figure_paired_difference_2d_*.pdf",
}
COLOUR_SEMANTICS = {
    "dominant": "#0072B2",
    "weighted": "#D55E00",
    "difference": "viridis colour-safe palette by alignment length",
}
ORIGINAL_REFERENCE_ROOT = CASE_ROOT / "original-paper"
ORIGINAL_REFERENCES = {
    "simulation-main.pdf": "references/SatuTe_manuscript_2024_08/Base_Nature_Methods/figures/figure_simulated_data_main.pdf",
    "simulation-misspecification.pdf": "references/SatuTe_manuscript_2024_08/Base_Nature_Methods/figures/figure_simulated_data_misspecification.pdf",
}


def link_data(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    relative_source = Path(os.path.relpath(source, destination.parent))
    if destination.is_symlink():
        if Path(os.readlink(destination)) == relative_source:
            return
        destination.unlink()
    elif destination.exists():
        raise RuntimeError(f"refusing to replace non-link data file: {destination}")
    destination.symlink_to(relative_source)


def prepare_case(case: Case) -> tuple[Path, Path, Path]:
    case_dir = CASE_ROOT / case.id
    source_dir = PROJECT / case.dataset
    summary = source_dir / "head_to_head_summary.tsv"
    detail = source_dir / "head_to_head_detail.tsv"
    for source in (summary, detail):
        if not source.is_file():
            raise FileNotFoundError(f"{case.id}: missing source data: {source}")
    link_data(summary, case_dir / "data" / "head_to_head_summary.tsv")
    link_data(detail, case_dir / "data" / "head_to_head_detail.tsv")
    (case_dir / "figures").mkdir(parents=True, exist_ok=True)
    return case_dir, summary, detail


def render_case(case: Case, case_dir: Path, summary: Path, detail: Path, bootstrap_reps: int) -> None:
    raw_dir = case_dir / "rendered"
    raw_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        (
            "Rscript",
            str(RENDERER),
            "--summary",
            str(summary),
            "--detail",
            str(detail),
            "--outdir",
            str(raw_dir),
            "--simulation-model",
            case.simulation_model,
            "--evaluation-model",
            case.evaluation_model,
            "--expected-reps",
            "1000",
            "--bootstrap-reps",
            str(bootstrap_reps),
            "--decision-rule",
            "taxon_bonferroni",
            "--scenario",
            "ml_tree",
        ),
        cwd=PROJECT,
        check=True,
    )
    for plot_kind, pattern in PLOT_FILES.items():
        matches = list(raw_dir.glob(pattern))
        if len(matches) != 1:
            raise RuntimeError(f"{case.id}: expected one {pattern}, found {len(matches)}")
        shutil.copy2(matches[0], case_dir / "figures" / f"{plot_kind}.pdf")
    paired_tables = list(raw_dir.glob("paired_formula_difference_*.tsv"))
    if len(paired_tables) != 1:
        raise RuntimeError(
            f"{case.id}: expected one paired_formula_difference_*.tsv, "
            f"found {len(paired_tables)}"
        )
    shutil.copy2(
        paired_tables[0],
        case_dir / "data" / "paired_formula_difference.tsv",
    )
    if case.id == "jc_JC":
        original_style = case_dir / "figures" / "original-style-baseline.pdf"
        subprocess.run(
            (
                "Rscript",
                os.fspath(BASELINE_RENDERER),
                "--summary",
                os.fspath(summary),
                "--output",
                os.fspath(original_style),
            ),
            cwd=PROJECT,
            check=True,
        )
        shutil.copy2(
            original_style,
            PROJECT / "manuscript" / "figures" / "figure_jc_baseline_power.pdf",
        )


def index_existing(case: Case, case_dir: Path) -> None:
    for plot_kind in PLOT_FILES:
        source = PUBLICATION_ROOT / f"{case.id}_{plot_kind}.pdf"
        if not source.is_file():
            raise FileNotFoundError(f"{case.id}: missing existing publication plot: {source}")
        shutil.copy2(source, case_dir / "figures" / f"{plot_kind}.pdf")
    if case.id == "jc_JC":
        shutil.copy2(
            PROJECT / "manuscript" / "figures" / "figure_jc_baseline_power.pdf",
            case_dir / "figures" / "original-style-baseline.pdf",
        )


def publish_case(case: Case, case_dir: Path) -> None:
    PUBLICATION_ROOT.mkdir(parents=True, exist_ok=True)
    for plot_kind in PLOT_FILES:
        source = case_dir / "figures" / f"{plot_kind}.pdf"
        if not source.is_file():
            raise FileNotFoundError(f"{case.id}: missing canonical plot: {source}")
        shutil.copy2(source, PUBLICATION_ROOT / f"{case.id}_{plot_kind}.pdf")

    provenance = {
        **asdict(case),
        "decision_rule": "taxon_bonferroni",
        "scenario": "ml_tree",
        "expected_replicates": 1000,
        "colour_semantics": COLOUR_SEMANTICS,
        "data_policy": "case-local symbolic links to immutable merged source tables",
    }
    if case.id == "jc_JC":
        provenance["additional_figure"] = (
            "figures/original-style-baseline.pdf: visual reconstruction of the "
            "original SatuTe two-panel JC figure with unrooted ape phylogenies "
            "rendered by ggtree"
        )
    (case_dir / "provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def index_original_paper() -> None:
    """Keep the archived figures beside explicit reconstruction data.

    The manuscript archive contains the final PDFs and methods text but no raw
    historical simulation tables. The current 1,000-replicate JC and
    misspecification tables are therefore linked as reconstruction data and are
    labelled as such rather than misrepresented as the original raw inputs.
    """

    for name, relative_source in ORIGINAL_REFERENCES.items():
        source = PROJECT / relative_source
        if not source.is_file():
            raise FileNotFoundError(f"missing original-paper reference: {source}")
        link_data(source, ORIGINAL_REFERENCE_ROOT / "figures" / name)

    reconstruction_data = {
        "jc_summary.tsv": PROJECT / JC_DATASET / "head_to_head_summary.tsv",
        "jc_detail.tsv": PROJECT / JC_DATASET / "head_to_head_detail.tsv",
        "misspecification_summary.tsv": PROJECT / MISSPEC_DATASET / "head_to_head_summary.tsv",
        "misspecification_detail.tsv": PROJECT / MISSPEC_DATASET / "head_to_head_detail.tsv",
        "archived_methods.tex": (
            PROJECT
            / "references/SatuTe_manuscript_2024_08/Base_Nature_Methods/methods_2024_08_01.tex"
        ),
    }
    for name, source in reconstruction_data.items():
        link_data(source, ORIGINAL_REFERENCE_ROOT / "data" / name)

    (ORIGINAL_REFERENCE_ROOT / "README.md").write_text(
        "# Original SatuTe figure references\n\n"
        "This directory preserves the two archived paper figures and links the "
        "current full-grid reconstruction tables used for comparison. The paper "
        "archive did not contain its historical raw simulation tables, so these "
        "links are explicitly reconstruction data, not relabelled original inputs.\n",
        encoding="utf-8",
    )


def render_integrated_main_figures() -> None:
    """Render and publish Figures 5 and 6 from stable case-local tables."""

    subprocess.run(
        (sys.executable, os.fspath(MISSPECIFICATION_RENDERER)),
        cwd=PROJECT,
        check=True,
    )
    subprocess.run(
        (sys.executable, os.fspath(CROSS_MODEL_RENDERER)),
        cwd=PROJECT,
        check=True,
    )
    for name in (
        "figure_model_misspecification.pdf",
        "figure_power_curve_summary.pdf",
    ):
        figure_id = (
            "model-misspecification"
            if name == "figure_model_misspecification.pdf"
            else "cross-model-power"
        )
        source = INTEGRATED_ROOT / figure_id / "figures" / name
        if not source.is_file():
            raise FileNotFoundError(f"integrated renderer did not create {source}")
        shutil.copy2(source, PROJECT / "manuscript" / "figures" / name)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", dest="cases", help="Build only this case; repeat as needed.")
    parser.add_argument("--mode", choices=("render", "index-existing"), default="render")
    parser.add_argument("--bootstrap-reps", type=int, default=5000)
    parser.add_argument("--plan", action="store_true")
    args = parser.parse_args()

    selected = [case for case in CASES if not args.cases or case.id in args.cases]
    unknown = sorted(set(args.cases or ()) - {case.id for case in CASES})
    if unknown:
        parser.error(f"unknown cases: {', '.join(unknown)}")
    if args.bootstrap_reps < 1000:
        parser.error("--bootstrap-reps must be at least 1000")

    print("Canonical publication power-figure plan:")
    for case in selected:
        print(f"  {case.id}: {case.simulation_model} -> {case.evaluation_model} [{case.dataset}]")
    if args.plan:
        return 0

    index_original_paper()
    for case in selected:
        print(f"\n[{case.id}] preparing case-local data and figures")
        case_dir, summary, detail = prepare_case(case)
        if args.mode == "render":
            render_case(case, case_dir, summary, detail, args.bootstrap_reps)
        else:
            index_existing(case, case_dir)
        publish_case(case, case_dir)

    selected_ids = {case.id for case in selected}
    integrated_dependencies = {
        "jc_JC",
        "gtr_GTR_PF06346",
        "gtr_GTR_PF06346_G4",
        "protein_LG_G4",
        "misspec_JC",
        "misspec_K2P",
        "misspec_F81",
    }
    if args.mode == "render" and integrated_dependencies <= selected_ids:
        print("\n[main figures] rendering contract-checked Figures 5 and 6")
        render_integrated_main_figures()

    print(f"\nPublished {len(selected)} complete power-figure cases from {CASE_ROOT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
