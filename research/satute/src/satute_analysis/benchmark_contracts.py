"""Stable schemas and decision-rule strategies for SatuTe benchmarks.

The benchmark deliberately separates the tree-analysis scenario from the
multiple-testing rule.  One native ``.sat.stat`` file supplies every decision
rule, so corrections never trigger duplicate simulations or IQ-TREE runs.
"""

from __future__ import annotations

from dataclasses import dataclass


SCHEMA_VERSION = 2
FORMULAS = ("dominant", "eigenvalue_weighted")

TRUE_FIXED = "true_tree_fixed_lengths"
TRUE_ML_LENGTHS = "true_topology_ml_lengths"
ML_TREE = "ml_tree"

FIG2_SCENARIOS = (TRUE_FIXED, TRUE_ML_LENGTHS, ML_TREE)
MISSPECIFICATION_SCENARIOS = (TRUE_FIXED, TRUE_ML_LENGTHS, ML_TREE)


@dataclass(frozen=True)
class DecisionRule:
    """Strategy describing how one native decision is summarized."""

    name: str
    decision_column: str
    p_value_column: str
    alpha_column: str


DECISION_RULES = (
    DecisionRule("unadjusted", "decision_unadjusted", "satP", "alpha"),
    DecisionRule(
        "taxon_bonferroni",
        "decision_taxon_bonf",
        "p_taxon_bonf",
        "alpha_taxon_bonf",
    ),
    DecisionRule("by_fdr", "decision_fdr", "fdr_by", "alpha"),
)
DECISION_RULE_BY_NAME = {rule.name: rule for rule in DECISION_RULES}


BASE_DETAIL_FIELDS = (
    "schema_version",
    "tree_case",
    "simulation_model",
    "evaluation_model",
    "fitted_evaluation_model",
    "nsites",
    "branch_length",
    "replicate",
    "seed",
    "simulator",
    "branch_length_source",
    "target_split",
    "scenario",
    "formula",
    "implementation",
    "target_found",
    "branch_id",
    "left_taxa",
    "right_taxa",
    "valid_sites",
    "skipped_sites",
    "branch_length_used",
    "alpha",
    "satC",
    "satVar",
    "satSE",
    "satZ",
    "satP",
    "decision_unadjusted",
    "alpha_taxon_bonf",
    "p_taxon_bonf",
    "decision_taxon_bonf",
    "fdr_by",
    "decision_fdr",
    "modes",
    "eigenvalues",
    "weights",
)

SUMMARY_FIELDS = (
    "schema_version",
    "tree_case",
    "simulation_model",
    "evaluation_model",
    "nsites",
    "branch_length",
    "scenario",
    "decision_rule",
    "formula",
    "evaluated",
    "informative",
    "missing_split",
    "fraction_informative",
)

FDR_AUDIT_FIELDS = (
    "schema_version",
    "tree_case",
    "simulation_model",
    "evaluation_model",
    "nsites",
    "branch_length",
    "replicate",
    "seed",
    "scenario",
    "formula",
    "family_size",
    "valid_p_values",
    "unadjusted_rejections",
    "by_fdr_rejections",
    "minimum_p_value",
    "minimum_fdr_by",
)

# One compact row per pooled native branch and formula.  This is deliberately
# separate from the target-branch detail table: a target curve is not enough to
# estimate a tree-wide FDR or false-discovery proportion.
BRANCH_AUDIT_FIELDS = (
    "schema_version",
    "tree_case",
    "simulation_model",
    "evaluation_model",
    "nsites",
    "branch_length",
    "replicate",
    "seed",
    "scenario",
    "formula",
    "branch_id",
    "split",
    "left_taxa",
    "right_taxa",
    "satP",
    "decision_unadjusted",
    "alpha_taxon_bonf",
    "p_taxon_bonf",
    "decision_taxon_bonf",
    "fdr_by",
    "decision_fdr",
    "is_target",
)

TIMING_FIELDS = (
    "schema_version",
    "tree_case",
    "simulation_model",
    "evaluation_model",
    "nsites",
    "branch_length",
    "replicate",
    "seed",
    "simulation_cache_hit",
    "simulation_seconds",
    "true_fixed_seconds",
    "true_ml_lengths_seconds",
    "ml_tree_seconds",
    "total_seconds",
)


def selected_scenarios(simulation_model: str, evaluation_model: str, scenario_set: str) -> list[str]:
    """Return native analysis scenarios, without duplicating decision rules."""

    if scenario_set in {"fig2", "all"}:
        return list(FIG2_SCENARIOS)
    if scenario_set == "misspecification":
        if simulation_model == evaluation_model:
            return [TRUE_FIXED]
        return list(MISSPECIFICATION_SCENARIOS)
    raise ValueError(f"Unknown scenario set: {scenario_set}")


def applicable_decision_rules(scenario: str) -> tuple[DecisionRule, ...]:
    """Return scientifically meaningful decisions for a native tree analysis.

    BY is a tree-wide branch scan and is available for every tree analysis.
    The paper's taxon-pair correction is retained for the inferred-tree search,
    where possible placements were searched.  It is not relabelled as FDR.
    """

    if scenario == ML_TREE:
        return DECISION_RULES
    return (DECISION_RULE_BY_NAME["unadjusted"], DECISION_RULE_BY_NAME["by_fdr"])


def expected_rows_for_task(simulation_model: str, evaluation_model: str, scenario_set: str) -> int:
    return len(selected_scenarios(simulation_model, evaluation_model, scenario_set)) * len(FORMULAS)
