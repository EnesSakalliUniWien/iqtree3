#!/usr/bin/env python3
"""Combine, validate, and summarize truth-labelled FDR calibration shards."""

from __future__ import annotations

import argparse
import csv
import itertools
import math
from collections import defaultdict
from pathlib import Path

from satute_analysis.fdr_calibration import (
    CALIBRATION_SUMMARY_FIELDS,
    REPLICATE_FDP_FIELDS,
    TRUTH_BRANCH_FIELDS,
    aggregate_replicate_fdp,
)


TIMING_FIELDS = (
    "schema_version", "design", "tree_case", "simulation_model", "nsites",
    "branch_length", "replicate", "seed", "truth", "simulation_seconds",
    "analysis_seconds", "total_seconds",
)
HELDOUT_FIELDS = (
    "schema_version", "design", "tree_case", "simulation_model",
    "evaluation_model", "nsites", "branch_length", "formula", "decision_rule",
    "training_nulls", "empirical_alpha", "empirical_threshold",
    "holdout_nulls", "holdout_rejections", "holdout_rejection_rate",
)


def read_table(path, fields):
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise ValueError(f"Unexpected schema in {path}: {reader.fieldnames}")
        return list(reader)


def write_table(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_manifest(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != ("field", "value"):
            raise ValueError(f"Unexpected manifest schema in {path}: {reader.fieldnames}")
        rows = list(reader)
    manifest = {row["field"]: row["value"] for row in rows}
    if len(manifest) != len(rows):
        raise ValueError(f"Duplicate field in {path}")
    return manifest


def csv_values(value, cast=str):
    return [cast(item.strip()) for item in value.split(",") if item.strip()]


def calibration_cell(row):
    return (
        row["design"], row["tree_case"], row["simulation_model"],
        row["evaluation_model"], str(int(row["nsites"])),
        f"{float(row['branch_length']):.12g}",
        row["formula"], row["decision_rule"],
    )


def validate_replicates(
    truth_rows,
    replicate_rows,
    expected_reps,
    training_reps,
    expected_cells,
):
    replicate_keys = set()
    by_cell = defaultdict(set)
    for row in replicate_rows:
        key = (*calibration_cell(row), row["replicate"])
        if key in replicate_keys:
            raise ValueError(f"Duplicate replicate FDP row: {key}")
        replicate_keys.add(key)
        by_cell[calibration_cell(row)].add(int(row["replicate"]))
    expected = set(range(1, expected_reps + 1))
    observed_cells = set(by_cell)
    if observed_cells != expected_cells:
        missing = sorted(expected_cells - observed_cells)[:10]
        extra = sorted(observed_cells - expected_cells)[:10]
        raise ValueError(f"Calibration-cell mismatch: missing={missing}, extra={extra}")
    for cell, replicates in by_cell.items():
        if replicates != expected:
            missing = sorted(expected - replicates)[:10]
            extra = sorted(replicates - expected)[:10]
            raise ValueError(f"Incomplete replicate cell {cell}: missing={missing}, extra={extra}")

    branch_counts = defaultdict(int)
    null_counts = defaultdict(int)
    truth_keys = set()
    for row in truth_rows:
        key = (*calibration_cell(row), row["replicate"])
        truth_key = (*key, row["split"])
        if truth_key in truth_keys:
            raise ValueError(f"Duplicate truth-labelled branch row: {truth_key}")
        truth_keys.add(truth_key)
        branch_counts[key] += 1
        null_counts[key] += int(row["truth"] == "null")
    for row in replicate_rows:
        key = (*calibration_cell(row), row["replicate"])
        if branch_counts[key] != int(row["family_size"]):
            raise ValueError(f"Truth-family size mismatch for {key}")
        if null_counts[key] != int(row["null_branches"]):
            raise ValueError(f"Truth null-count mismatch for {key}")
        if row["design"] == "exact_null":
            expected_sample = "training" if int(row["replicate"]) <= training_reps else "holdout"
            if row["sample"] != expected_sample or int(row["null_branches"]) != 1:
                raise ValueError(f"Exact-null sample/truth mismatch for {key}: {row}")
        elif row["design"] == "mixed_null":
            expected_nulls = 1 if int(row["replicate"]) % 2 else 0
            if row["sample"] != "operating" or int(row["null_branches"]) != expected_nulls:
                raise ValueError(f"Mixed-null prespecification mismatch for {key}: {row}")


def heldout_calibration(truth_rows, alpha):
    training = defaultdict(list)
    holdout = defaultdict(list)
    for row in truth_rows:
        if row["design"] != "exact_null" or row["truth"] != "null":
            continue
        try:
            p_value = float(row["p_value"])
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Null branch has invalid p-value: {row}") from exc
        target = training if row["sample"] == "training" else holdout
        target[calibration_cell(row)].append(p_value)

    rows = []
    for cell in sorted(training):
        values = sorted(training[cell])
        heldout = holdout.get(cell, [])
        if not values or not heldout:
            raise ValueError(f"Missing training or holdout nulls for {cell}")
        rank = max(1, min(len(values), math.floor(alpha * (len(values) + 1))))
        threshold = values[rank - 1]
        rejected = sum(value <= threshold for value in heldout)
        design, tree_case, simulation_model, evaluation_model, nsites, branch_length, formula, rule = cell
        rows.append(
            {
                "schema_version": 2,
                "design": design,
                "tree_case": tree_case,
                "simulation_model": simulation_model,
                "evaluation_model": evaluation_model,
                "nsites": nsites,
                "branch_length": branch_length,
                "formula": formula,
                "decision_rule": rule,
                "training_nulls": len(values),
                "empirical_alpha": alpha,
                "empirical_threshold": threshold,
                "holdout_nulls": len(heldout),
                "holdout_rejections": rejected,
                "holdout_rejection_rate": rejected / len(heldout),
            }
        )
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--expected-shards", type=int, default=1000)
    parser.add_argument("--expected-reps", type=int, required=True)
    parser.add_argument("--training-reps", type=int, default=0)
    parser.add_argument("--alpha", type=float, default=0.05)
    args = parser.parse_args()

    root = Path(args.root)
    outdir = Path(args.outdir)
    if not 0.0 < args.alpha < 1.0:
        raise SystemExit("--alpha must be strictly between zero and one")
    if args.expected_reps < 1 or args.expected_shards < 1:
        raise SystemExit("--expected-reps and --expected-shards must be positive")
    manifest = read_manifest(root / "run_manifest.tsv")
    required_manifest = {
        "schema_version",
        "design",
        "replicates",
        "training_replicates",
        "shard_count",
        "models",
        "resolved_model_pairs",
        "site_lengths",
        "branch_lengths",
        "tree_cases",
    }
    missing_manifest = sorted(required_manifest - set(manifest))
    if missing_manifest:
        raise SystemExit(f"Manifest is missing required fields: {missing_manifest}")
    if int(manifest["schema_version"]) != 2:
        raise SystemExit(f"Unsupported manifest schema: {manifest['schema_version']}")
    if int(manifest["replicates"]) != args.expected_reps:
        raise SystemExit("Manifest replicate count does not match --expected-reps")
    if int(manifest["training_replicates"]) != args.training_reps:
        raise SystemExit("Manifest training count does not match --training-reps")
    if int(manifest["shard_count"]) != args.expected_shards:
        raise SystemExit("Manifest shard count does not match --expected-shards")
    if manifest["design"] == "exact_null":
        if not 0 < args.training_reps < args.expected_reps:
            raise SystemExit("Exact-null postprocessing requires a non-empty training/holdout split")
    elif manifest["design"] == "mixed_null":
        if args.training_reps != 0 or args.expected_reps % 2:
            raise SystemExit("Mixed-null postprocessing requires zero training reps and even reps")
    else:
        raise SystemExit(f"Unsupported calibration design: {manifest['design']}")
    models = csv_values(manifest["models"])
    sites = csv_values(manifest["site_lengths"], int)
    branches = csv_values(manifest["branch_lengths"], float)
    trees = csv_values(manifest["tree_cases"])
    resolved_pairs = manifest["resolved_model_pairs"].split(";")
    if len(resolved_pairs) != len(models) or any("=>" not in pair for pair in resolved_pairs):
        raise SystemExit("Manifest resolved-model pairs do not match the model list")
    expected_cells = {
        (
            manifest["design"],
            tree,
            model,
            model,
            str(site),
            f"{branch:.12g}",
            formula,
            rule,
        )
        for tree, model, site, branch, formula, rule in itertools.product(
            trees,
            models,
            sites,
            branches,
            ("dominant", "eigenvalue_weighted"),
            ("unadjusted", "taxon_bonferroni", "by_fdr"),
        )
    }
    outdir.mkdir(parents=True, exist_ok=True)
    shards = sorted(root.glob("shards/shard-*"))
    if len(shards) != args.expected_shards:
        raise SystemExit(f"Found {len(shards)} shard directories; expected {args.expected_shards}")

    truth_rows = []
    replicate_rows = []
    timing_rows = []
    for shard in shards:
        truth_rows.extend(read_table(shard / "truth_branch_audit.tsv", TRUTH_BRANCH_FIELDS))
        replicate_rows.extend(read_table(shard / "replicate_fdp.tsv", REPLICATE_FDP_FIELDS))
        timing_rows.extend(read_table(shard / "timing.tsv", TIMING_FIELDS))
    validate_replicates(
        truth_rows,
        replicate_rows,
        args.expected_reps,
        args.training_reps,
        expected_cells,
    )
    summaries = aggregate_replicate_fdp(replicate_rows)
    holdout = heldout_calibration(truth_rows, args.alpha) if args.training_reps else []

    write_table(outdir / "truth_branch_audit.tsv", TRUTH_BRANCH_FIELDS, truth_rows)
    write_table(outdir / "replicate_fdp.tsv", REPLICATE_FDP_FIELDS, replicate_rows)
    write_table(outdir / "calibration_summary.tsv", CALIBRATION_SUMMARY_FIELDS, summaries)
    write_table(outdir / "timing.tsv", TIMING_FIELDS, timing_rows)
    if holdout:
        write_table(outdir / "heldout_null_calibration.tsv", HELDOUT_FIELDS, holdout)
    print(f"Shards: {len(shards)}")
    print(f"Truth rows: {len(truth_rows)}")
    print(f"Replicate FDP rows: {len(replicate_rows)}")
    print(f"Summary rows: {len(summaries)}")
    print(f"Held-out rows: {len(holdout)}")


if __name__ == "__main__":
    main()
