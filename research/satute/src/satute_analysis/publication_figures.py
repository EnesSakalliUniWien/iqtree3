"""Scientific contracts and stable data paths for publication figures.

Renderers consume this module instead of reconstructing scenario/rule choices
from whichever columns happen to be present in a detail table.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .benchmark_contracts import ML_TREE, TRUE_FIXED, TRUE_ML_LENGTHS


PROJECT_ROOT = Path(__file__).resolve().parents[2]
PUBLICATION_CASE_ROOT = PROJECT_ROOT / "artifacts" / "publication" / "power"
PUBLICATION_FIGURE_ROOT = PROJECT_ROOT / "artifacts" / "publication" / "figures"


@dataclass(frozen=True)
class FigureSeriesContract:
    id: str
    label: str
    scenario: str
    decision_rule: str
    reference_model: bool = False


@dataclass(frozen=True)
class MisspecificationFigureContract:
    id: str
    series: tuple[FigureSeriesContract, ...]
    evaluation_models: tuple[str, ...]
    tree_cases: tuple[str, ...]


@dataclass(frozen=True)
class CrossModelFigureContract:
    id: str
    scenario: str
    decision_rule: str
    site_lengths: tuple[int, ...]
    case_ids: tuple[str, ...]


MISSPECIFICATION_CONTRACT = MisspecificationFigureContract(
    id="model-misspecification",
    series=(
        FigureSeriesContract(
            "correct_reference",
            "Correct GTR+F, true tree",
            TRUE_FIXED,
            "unadjusted",
            reference_model=True,
        ),
        FigureSeriesContract(
            "misspecified_true_tree",
            "Misspecified model, true tree",
            TRUE_FIXED,
            "unadjusted",
        ),
        FigureSeriesContract(
            "misspecified_ml_lengths",
            "Misspecified model, true topology and ML lengths",
            TRUE_ML_LENGTHS,
            "unadjusted",
        ),
        FigureSeriesContract(
            "misspecified_ml_tree_bonferroni",
            "Misspecified model, ML tree and Bonferroni",
            ML_TREE,
            "taxon_bonferroni",
        ),
    ),
    evaluation_models=("JC", "K2P", "F81"),
    tree_cases=("five_external", "sixteen_internal"),
)


CROSS_MODEL_CONTRACT = CrossModelFigureContract(
    id="cross-model-power",
    scenario=ML_TREE,
    decision_rule="taxon_bonferroni",
    site_lengths=(100, 1000),
    case_ids=("jc_JC", "gtr_GTR_PF06346_G4", "protein_LG_G4"),
)


def stable_case_summary(case_id: str) -> Path:
    """Return the stable case-local summary path used by publication renderers."""

    return PUBLICATION_CASE_ROOT / case_id / "data" / "head_to_head_summary.tsv"


def stable_figure_root(figure_id: str) -> Path:
    """Return the self-contained publication directory for one figure."""

    return PUBLICATION_FIGURE_ROOT / figure_id
