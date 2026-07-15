#!/usr/bin/env python3

from satute_analysis.fdr_calibration import aggregate_replicate_fdp, truth_labelled_family


def branch(formula, branch_id, split, raw, bonf, by):
    return {
        "formula": formula,
        "branch_id": branch_id,
        "split": split,
        "satP": raw,
        "decision_unadjusted": "informative" if raw <= 0.05 else "saturated",
        "p_taxon_bonf": bonf,
        "decision_taxon_bonf": "informative" if bonf <= 0.05 else "saturated",
        "fdr_by": by,
        "decision_fdr": "informative" if by <= 0.05 else "saturated",
    }


def main():
    rows = []
    for formula in ("dominant", "eigenvalue_weighted"):
        rows.extend(
            (
                branch(formula, "1", "B", 0.01, 0.06, 0.08),
                branch(formula, "2", "A1", 0.001, 0.004, 0.006),
                branch(formula, "3", "A2", 0.8, 1.0, 1.0),
            )
        )
    base = {
        "design": "exact_null",
        "sample": "holdout",
        "tree_case": "five_external",
        "simulation_model": "GTR",
        "evaluation_model": "GTR",
        "nsites": 100,
        "branch_length": 8,
        "replicate": 1,
        "seed": 17,
    }
    truth, replicate = truth_labelled_family(rows, [("B",)], base)
    if len(truth) != 18 or len(replicate) != 6:
        raise SystemExit(f"Unexpected calibration row counts: {len(truth)}, {len(replicate)}")
    dominant_raw = next(
        row for row in replicate
        if row["formula"] == "dominant" and row["decision_rule"] == "unadjusted"
    )
    if (dominant_raw["null_branches"], dominant_raw["false_discoveries"], dominant_raw["true_discoveries"]) != (1, 1, 1):
        raise SystemExit(f"Incorrect truth-labelled counts: {dominant_raw}")
    dominant_by = next(
        row for row in replicate
        if row["formula"] == "dominant" and row["decision_rule"] == "by_fdr"
    )
    if (dominant_by["false_discoveries"], dominant_by["true_discoveries"]) != (0, 1):
        raise SystemExit(f"Incorrect BY counts: {dominant_by}")
    summary = aggregate_replicate_fdp(replicate)
    raw_summary = next(
        row for row in summary
        if row["formula"] == "dominant" and row["decision_rule"] == "unadjusted"
    )
    if (
        raw_summary["mean_fdp"] != 0.5
        or raw_summary["pooled_false_discovery_fraction"] != 0.5
        or raw_summary["null_rejection_rate"] != 1.0
    ):
        raise SystemExit(f"Incorrect aggregate calibration summary: {raw_summary}")
    try:
        truth_labelled_family(rows + [rows[0]], [("B",)], base)
    except ValueError as exc:
        if "Duplicate normalized splits" not in str(exc):
            raise
    else:
        raise SystemExit("Duplicate split must fail the truth-map gate")
    try:
        truth_labelled_family(rows, [("missing",)], base)
    except ValueError as exc:
        if "null splits absent" not in str(exc):
            raise
    else:
        raise SystemExit("An absent generative null split must fail the truth-map gate")
    invalid_rows = [dict(row) for row in rows]
    invalid_rows[0]["satP"] = ""
    try:
        truth_labelled_family(invalid_rows, [("B",)], base)
    except ValueError as exc:
        if "Invalid satP" not in str(exc):
            raise
    else:
        raise SystemExit("A missing calibration p-value must fail the audit")
    print("Truth-labelled FDR calibration checks passed")


if __name__ == "__main__":
    main()
