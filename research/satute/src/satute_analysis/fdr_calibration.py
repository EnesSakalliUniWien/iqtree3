"""Truth-labelled branch-family summaries for SatuTe FDR calibration."""

from __future__ import annotations

from collections import defaultdict

from .benchmark_contracts import SCHEMA_VERSION
from .native_results import normalized_split


CALIBRATION_RULES = (
    ("unadjusted", "decision_unadjusted", "satP"),
    ("taxon_bonferroni", "decision_taxon_bonf", "p_taxon_bonf"),
    ("by_fdr", "decision_fdr", "fdr_by"),
)

TRUTH_BRANCH_FIELDS = (
    "schema_version",
    "design",
    "sample",
    "tree_case",
    "simulation_model",
    "evaluation_model",
    "nsites",
    "branch_length",
    "replicate",
    "seed",
    "formula",
    "branch_id",
    "split",
    "truth",
    "decision_rule",
    "p_value",
    "decision",
    "discovery",
    "false_discovery",
    "true_discovery",
)

REPLICATE_FDP_FIELDS = (
    "schema_version",
    "design",
    "sample",
    "tree_case",
    "simulation_model",
    "evaluation_model",
    "nsites",
    "branch_length",
    "replicate",
    "seed",
    "formula",
    "decision_rule",
    "family_size",
    "null_branches",
    "alternative_branches",
    "discoveries",
    "false_discoveries",
    "true_discoveries",
    "fdp",
    "power",
)

CALIBRATION_SUMMARY_FIELDS = (
    "schema_version",
    "design",
    "sample",
    "tree_case",
    "simulation_model",
    "evaluation_model",
    "nsites",
    "branch_length",
    "formula",
    "decision_rule",
    "replicates",
    "mean_fdp",
    "pooled_false_discovery_fraction",
    "null_rejection_rate",
    "power",
    "mean_discoveries",
)


def _probability(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return ""
    return number if 0.0 <= number <= 1.0 else ""


def truth_labelled_family(branch_rows, null_splits, base):
    """Label every pooled branch by split and emit rule-specific audit rows.

    ``null_splits`` is independent design metadata. No label is inferred from a
    p-value, adjusted p-value, or decision.
    """

    normalized_nulls = {normalized_split(split) for split in null_splits}
    by_formula = defaultdict(list)
    for row in branch_rows:
        split = normalized_split(row.get("split", ""))
        if not split:
            raise ValueError(f"Branch has no normalized split: {row}")
        by_formula[row["formula"]].append((split, row))

    truth_rows = []
    replicate_rows = []
    for formula, entries in sorted(by_formula.items()):
        splits = [split for split, _row in entries]
        if len(splits) != len(set(splits)):
            raise ValueError(f"Duplicate normalized splits in {formula} branch family")
        missing_nulls = normalized_nulls - set(splits)
        if missing_nulls:
            raise ValueError(
                f"Truth-labelled null splits absent from {formula} branch family: "
                f"{sorted(missing_nulls)}"
            )
        for rule_name, decision_column, p_column in CALIBRATION_RULES:
            labelled = []
            for split, row in entries:
                truth = "null" if split in normalized_nulls else "alternative"
                decision = row.get(decision_column, "")
                if decision not in {"informative", "saturated"}:
                    raise ValueError(
                        f"Invalid {decision_column}={decision!r} for {formula} split {split}"
                    )
                discovery = int(decision == "informative")
                false_discovery = int(discovery and truth == "null")
                true_discovery = int(discovery and truth == "alternative")
                p_value = _probability(row.get(p_column, ""))
                if p_value == "":
                    raise ValueError(
                        f"Invalid {p_column}={row.get(p_column)!r} for "
                        f"{formula} split {split}"
                    )
                output = {
                    "schema_version": SCHEMA_VERSION,
                    **base,
                    "formula": formula,
                    "branch_id": row.get("branch_id", ""),
                    "split": ",".join(split),
                    "truth": truth,
                    "decision_rule": rule_name,
                    "p_value": p_value,
                    "decision": decision,
                    "discovery": discovery,
                    "false_discovery": false_discovery,
                    "true_discovery": true_discovery,
                }
                truth_rows.append(output)
                labelled.append(output)

            null_count = sum(row["truth"] == "null" for row in labelled)
            alternative_count = sum(row["truth"] == "alternative" for row in labelled)
            discoveries = sum(row["discovery"] for row in labelled)
            false_discoveries = sum(row["false_discovery"] for row in labelled)
            true_discoveries = sum(row["true_discovery"] for row in labelled)
            replicate_rows.append(
                {
                    "schema_version": SCHEMA_VERSION,
                    **base,
                    "formula": formula,
                    "decision_rule": rule_name,
                    "family_size": len(labelled),
                    "null_branches": null_count,
                    "alternative_branches": alternative_count,
                    "discoveries": discoveries,
                    "false_discoveries": false_discoveries,
                    "true_discoveries": true_discoveries,
                    "fdp": false_discoveries / max(1, discoveries),
                    "power": true_discoveries / alternative_count if alternative_count else "",
                }
            )
    return truth_rows, replicate_rows


def aggregate_replicate_fdp(rows):
    grouped = defaultdict(list)
    key_fields = (
        "design",
        "sample",
        "tree_case",
        "simulation_model",
        "evaluation_model",
        "nsites",
        "branch_length",
        "formula",
        "decision_rule",
    )
    for row in rows:
        grouped[tuple(row[field] for field in key_fields)].append(row)

    summaries = []
    for key, values in sorted(grouped.items()):
        replicates = len(values)
        false_discoveries = sum(int(row["false_discoveries"]) for row in values)
        discoveries = sum(int(row["discoveries"]) for row in values)
        null_branches = sum(int(row["null_branches"]) for row in values)
        true_discoveries = sum(int(row["true_discoveries"]) for row in values)
        alternative_branches = sum(int(row["alternative_branches"]) for row in values)
        summaries.append(
            {
                "schema_version": SCHEMA_VERSION,
                **dict(zip(key_fields, key)),
                "replicates": replicates,
                "mean_fdp": sum(float(row["fdp"]) for row in values) / replicates,
                "pooled_false_discovery_fraction": false_discoveries / max(1, discoveries),
                "null_rejection_rate": false_discoveries / null_branches if null_branches else "",
                "power": true_discoveries / alternative_branches if alternative_branches else "",
                "mean_discoveries": discoveries / replicates,
            }
        )
    return summaries
