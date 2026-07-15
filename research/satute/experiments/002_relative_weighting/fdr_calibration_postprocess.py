#!/usr/bin/env python3
"""Stream, validate, and summarize truth-labelled FDR calibration shards."""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import os
from collections import defaultdict
from pathlib import Path

from satute_analysis.benchmark_contracts import SCHEMA_VERSION
from satute_analysis.fdr_calibration import (
    CALIBRATION_SUMMARY_FIELDS,
    REPLICATE_FDP_FIELDS,
    TRUTH_BRANCH_FIELDS,
)


TIMING_FIELDS = (
    "schema_version",
    "design",
    "tree_case",
    "simulation_model",
    "nsites",
    "branch_length",
    "replicate",
    "seed",
    "truth",
    "simulation_seconds",
    "analysis_seconds",
    "total_seconds",
)
HELDOUT_FIELDS = (
    "schema_version",
    "design",
    "tree_case",
    "simulation_model",
    "evaluation_model",
    "nsites",
    "branch_length",
    "formula",
    "decision_rule",
    "training_nulls",
    "empirical_alpha",
    "empirical_threshold",
    "holdout_nulls",
    "holdout_rejections",
    "holdout_rejection_rate",
)
FORMULAS = ("dominant", "eigenvalue_weighted")
RULES = ("unadjusted", "taxon_bonferroni", "by_fdr")
SUMMARY_KEY_FIELDS = (
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


def iter_table(path, fields):
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise ValueError(f"Unexpected schema in {path}: {reader.fieldnames}")
        yield from reader


def validate_header(path, fields):
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
    if tuple(header or ()) != tuple(fields):
        raise ValueError(f"Unexpected schema in {path}: {header}")


def open_writer(path, fields):
    handle = path.open("w", encoding="utf-8", newline="")
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    return handle, writer


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


def canonical_int(value):
    return str(int(value))


def canonical_float(value):
    return f"{float(value):.12g}"


def probability(value, label):
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid {label}: {value!r}") from exc
    if not 0.0 <= result <= 1.0:
        raise ValueError(f"Invalid {label}: {value!r}")
    return result


def calibration_cell(row):
    return (
        row["design"],
        row["tree_case"],
        row["simulation_model"],
        row["evaluation_model"],
        canonical_int(row["nsites"]),
        canonical_float(row["branch_length"]),
        row["formula"],
        row["decision_rule"],
    )


def replicate_key(row):
    return (*calibration_cell(row), canonical_int(row["replicate"]))


def expected_sample_and_nulls(design, replicate, training_reps):
    if design == "exact_null":
        return ("training" if replicate <= training_reps else "holdout"), 1
    if design == "mixed_null":
        return "operating", (1 if replicate % 2 else 0)
    raise ValueError(f"Unsupported calibration design: {design}")


def validate_replicate_row(row, training_reps):
    if int(row["schema_version"]) != SCHEMA_VERSION:
        raise ValueError(f"Unsupported replicate schema: {row['schema_version']}")
    if row["formula"] not in FORMULAS or row["decision_rule"] not in RULES:
        raise ValueError(f"Unexpected formula/rule: {row}")
    replicate = int(row["replicate"])
    expected_sample, expected_nulls = expected_sample_and_nulls(
        row["design"], replicate, training_reps
    )
    family_size = int(row["family_size"])
    nulls = int(row["null_branches"])
    alternatives = int(row["alternative_branches"])
    discoveries = int(row["discoveries"])
    false_discoveries = int(row["false_discoveries"])
    true_discoveries = int(row["true_discoveries"])
    if row["sample"] != expected_sample or nulls != expected_nulls:
        raise ValueError(f"Sample/truth prespecification mismatch: {row}")
    if family_size <= 0 or nulls + alternatives != family_size:
        raise ValueError(f"Invalid branch-family denominator: {row}")
    if false_discoveries + true_discoveries != discoveries:
        raise ValueError(f"Discovery decomposition mismatch: {row}")
    if false_discoveries > nulls or true_discoveries > alternatives:
        raise ValueError(f"Discovery count exceeds its truth denominator: {row}")
    expected_fdp = false_discoveries / max(1, discoveries)
    if not math.isclose(float(row["fdp"]), expected_fdp, abs_tol=1e-12):
        raise ValueError(f"Replicate FDP mismatch: {row}")
    if alternatives:
        expected_power = true_discoveries / alternatives
        if not math.isclose(float(row["power"]), expected_power, abs_tol=1e-12):
            raise ValueError(f"Replicate power mismatch: {row}")
    elif row["power"] != "":
        raise ValueError(f"Power must be empty without alternatives: {row}")


def validate_truth_row(row, replicate_row):
    if int(row["schema_version"]) != SCHEMA_VERSION:
        raise ValueError(f"Unsupported truth schema: {row['schema_version']}")
    for field in (
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
    ):
        left = row[field]
        right = replicate_row[field]
        if field == "branch_length":
            equal = canonical_float(left) == canonical_float(right)
        elif field in {"nsites", "replicate", "seed"}:
            equal = canonical_int(left) == canonical_int(right)
        else:
            equal = left == right
        if not equal:
            raise ValueError(f"Truth/replicate field mismatch for {field}: {row}")
    if row["truth"] not in {"null", "alternative"}:
        raise ValueError(f"Invalid truth label: {row}")
    if row["decision"] not in {"informative", "saturated"}:
        raise ValueError(f"Invalid decision label: {row}")
    probability(row["p_value"], "truth-row p-value")
    discovery = int(row["discovery"])
    false_discovery = int(row["false_discovery"])
    true_discovery = int(row["true_discovery"])
    if discovery not in {0, 1} or false_discovery not in {0, 1} or true_discovery not in {0, 1}:
        raise ValueError(f"Truth decision flags must be binary: {row}")
    if discovery != int(row["decision"] == "informative"):
        raise ValueError(f"Truth decision flag mismatch: {row}")
    if false_discovery != int(discovery and row["truth"] == "null"):
        raise ValueError(f"False-discovery truth mismatch: {row}")
    if true_discovery != int(discovery and row["truth"] == "alternative"):
        raise ValueError(f"True-discovery truth mismatch: {row}")


def update_aggregate(aggregates, row):
    key = (
        row["design"],
        row["sample"],
        row["tree_case"],
        row["simulation_model"],
        row["evaluation_model"],
        canonical_int(row["nsites"]),
        canonical_float(row["branch_length"]),
        row["formula"],
        row["decision_rule"],
    )
    state = aggregates[key]
    state["replicates"] += 1
    state["sum_fdp"] += float(row["fdp"])
    for field in (
        "false_discoveries",
        "discoveries",
        "null_branches",
        "true_discoveries",
        "alternative_branches",
    ):
        state[field] += int(row[field])


def aggregate_rows(aggregates):
    rows = []
    for key, state in sorted(aggregates.items()):
        replicates = state["replicates"]
        discoveries = state["discoveries"]
        nulls = state["null_branches"]
        alternatives = state["alternative_branches"]
        rows.append(
            {
                "schema_version": SCHEMA_VERSION,
                **dict(zip(SUMMARY_KEY_FIELDS, key)),
                "replicates": replicates,
                "mean_fdp": state["sum_fdp"] / replicates,
                "pooled_false_discovery_fraction": state["false_discoveries"]
                / max(1, discoveries),
                "null_rejection_rate": state["false_discoveries"] / nulls if nulls else "",
                "power": state["true_discoveries"] / alternatives if alternatives else "",
                "mean_discoveries": discoveries / replicates,
            }
        )
    return rows


def heldout_calibration(training, holdout, alpha):
    rows = []
    if set(training) != set(holdout):
        raise ValueError("Training and held-out null calibration cells differ")
    for cell in sorted(training):
        values = sorted(training[cell])
        heldout = holdout[cell]
        if not values or not heldout:
            raise ValueError(f"Missing training or held-out nulls for {cell}")
        rank = max(1, min(len(values), math.floor(alpha * (len(values) + 1))))
        threshold = values[rank - 1]
        rejected = sum(value <= threshold for value in heldout)
        (
            design,
            tree_case,
            simulation_model,
            evaluation_model,
            nsites,
            branch_length,
            formula,
            rule,
        ) = cell
        rows.append(
            {
                "schema_version": SCHEMA_VERSION,
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


def expected_cells_from_manifest(manifest):
    models = csv_values(manifest["models"])
    sites = csv_values(manifest["site_lengths"], int)
    branches = csv_values(manifest["branch_lengths"], float)
    trees = csv_values(manifest["tree_cases"])
    return {
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
            trees, models, sites, branches, FORMULAS, RULES
        )
    }


def validate_manifest(manifest, args):
    required = {
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
    missing = sorted(required - set(manifest))
    if missing:
        raise ValueError(f"Manifest is missing required fields: {missing}")
    if int(manifest["schema_version"]) != SCHEMA_VERSION:
        raise ValueError(f"Unsupported manifest schema: {manifest['schema_version']}")
    if int(manifest["replicates"]) != args.expected_reps:
        raise ValueError("Manifest replicate count does not match --expected-reps")
    if int(manifest["training_replicates"]) != args.training_reps:
        raise ValueError("Manifest training count does not match --training-reps")
    if int(manifest["shard_count"]) != args.expected_shards:
        raise ValueError("Manifest shard count does not match --expected-shards")
    models = csv_values(manifest["models"])
    resolved_pairs = manifest["resolved_model_pairs"].split(";")
    if len(resolved_pairs) != len(models) or any("=>" not in pair for pair in resolved_pairs):
        raise ValueError("Manifest resolved-model pairs do not match the model list")
    if manifest["design"] == "exact_null":
        if not 0 < args.training_reps < args.expected_reps:
            raise ValueError("Exact-null postprocessing requires a training/holdout split")
    elif manifest["design"] == "mixed_null":
        if args.training_reps != 0 or args.expected_reps % 2:
            raise ValueError("Mixed-null postprocessing requires zero training reps and even reps")
    else:
        raise ValueError(f"Unsupported calibration design: {manifest['design']}")


def validate_global_completeness(by_cell, expected_cells, expected_reps):
    if set(by_cell) != expected_cells:
        missing = sorted(expected_cells - set(by_cell))[:10]
        extra = sorted(set(by_cell) - expected_cells)[:10]
        raise ValueError(f"Calibration-cell mismatch: missing={missing}, extra={extra}")
    expected = set(range(1, expected_reps + 1))
    for cell, replicates in by_cell.items():
        if replicates != expected:
            missing = sorted(expected - replicates)[:10]
            extra = sorted(replicates - expected)[:10]
            raise ValueError(f"Incomplete replicate cell {cell}: missing={missing}, extra={extra}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--expected-shards", type=int, default=1000)
    parser.add_argument("--expected-reps", type=int, required=True)
    parser.add_argument("--training-reps", type=int, default=0)
    parser.add_argument("--alpha", type=float, default=0.05)
    args = parser.parse_args()

    if not 0.0 < args.alpha < 1.0:
        raise SystemExit("--alpha must be strictly between zero and one")
    if args.expected_reps < 1 or args.expected_shards < 1:
        raise SystemExit("--expected-reps and --expected-shards must be positive")
    root = Path(args.root)
    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit(f"Refusing to overwrite non-empty output directory: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    manifest = read_manifest(root / "run_manifest.tsv")
    validate_manifest(manifest, args)
    expected_cells = expected_cells_from_manifest(manifest)
    shards = sorted(root.glob("shards/shard-*"))
    if len(shards) != args.expected_shards:
        raise SystemExit(f"Found {len(shards)} shard directories; expected {args.expected_shards}")

    output_names = {
        "truth": "truth_branch_audit.tsv",
        "replicate": "replicate_fdp.tsv",
        "summary": "calibration_summary.tsv",
        "timing": "timing.tsv",
        "heldout": "heldout_null_calibration.tsv",
    }
    temporary = {
        key: outdir / f".{name}.tmp.{os.getpid()}" for key, name in output_names.items()
    }
    truth_handle = replicate_handle = timing_handle = None
    truth_rows = replicate_rows = timing_rows = 0
    aggregates = defaultdict(
        lambda: {
            "replicates": 0,
            "sum_fdp": 0.0,
            "false_discoveries": 0,
            "discoveries": 0,
            "null_branches": 0,
            "true_discoveries": 0,
            "alternative_branches": 0,
        }
    )
    by_cell = defaultdict(set)
    timing_by_cell = defaultdict(set)
    training_nulls = defaultdict(list)
    holdout_nulls = defaultdict(list)
    published = []
    try:
        truth_handle, truth_writer = open_writer(temporary["truth"], TRUTH_BRANCH_FIELDS)
        replicate_handle, replicate_writer = open_writer(
            temporary["replicate"], REPLICATE_FDP_FIELDS
        )
        timing_handle, timing_writer = open_writer(temporary["timing"], TIMING_FIELDS)

        for shard_index, shard in enumerate(shards, start=1):
            validate_header(shard / "calibration_summary.tsv", CALIBRATION_SUMMARY_FIELDS)
            local_replicates = {}
            local_counts = defaultdict(
                lambda: {
                    "branches": 0,
                    "nulls": 0,
                    "discoveries": 0,
                    "false_discoveries": 0,
                    "true_discoveries": 0,
                }
            )
            local_splits = defaultdict(set)
            for row in iter_table(shard / "replicate_fdp.tsv", REPLICATE_FDP_FIELDS):
                validate_replicate_row(row, args.training_reps)
                key = replicate_key(row)
                if key in local_replicates or int(row["replicate"]) in by_cell[calibration_cell(row)]:
                    raise ValueError(f"Duplicate replicate FDP row: {key}")
                local_replicates[key] = row
                by_cell[calibration_cell(row)].add(int(row["replicate"]))
                update_aggregate(aggregates, row)
                replicate_writer.writerow(row)
                replicate_rows += 1

            for row in iter_table(shard / "truth_branch_audit.tsv", TRUTH_BRANCH_FIELDS):
                key = replicate_key(row)
                replicate_row = local_replicates.get(key)
                if replicate_row is None:
                    raise ValueError(f"Truth row has no shard-local replicate row: {key}")
                validate_truth_row(row, replicate_row)
                split = row["split"]
                if split in local_splits[key]:
                    raise ValueError(f"Duplicate truth-labelled split for {key}: {split}")
                local_splits[key].add(split)
                counts = local_counts[key]
                counts["branches"] += 1
                counts["nulls"] += int(row["truth"] == "null")
                counts["discoveries"] += int(row["discovery"])
                counts["false_discoveries"] += int(row["false_discovery"])
                counts["true_discoveries"] += int(row["true_discovery"])
                if row["design"] == "exact_null" and row["truth"] == "null":
                    destination = (
                        training_nulls if row["sample"] == "training" else holdout_nulls
                    )
                    destination[calibration_cell(row)].append(float(row["p_value"]))
                truth_writer.writerow(row)
                truth_rows += 1

            for key, replicate_row in local_replicates.items():
                counts = local_counts[key]
                comparisons = {
                    "branches": "family_size",
                    "nulls": "null_branches",
                    "discoveries": "discoveries",
                    "false_discoveries": "false_discoveries",
                    "true_discoveries": "true_discoveries",
                }
                for count_name, row_name in comparisons.items():
                    if counts[count_name] != int(replicate_row[row_name]):
                        raise ValueError(
                            f"Truth-family {count_name} mismatch for {key}: "
                            f"{counts[count_name]} != {replicate_row[row_name]}"
                        )

            for row in iter_table(shard / "timing.tsv", TIMING_FIELDS):
                if int(row["schema_version"]) != SCHEMA_VERSION:
                    raise ValueError(f"Unsupported timing schema in {shard}")
                replicate = int(row["replicate"])
                expected_sample, expected_nulls = expected_sample_and_nulls(
                    row["design"], replicate, args.training_reps
                )
                expected_truth = "null" if expected_nulls else "alternative"
                if row["truth"] != expected_truth:
                    raise ValueError(f"Timing truth mismatch: {row}")
                timing_cell = (
                    row["design"],
                    row["tree_case"],
                    row["simulation_model"],
                    canonical_int(row["nsites"]),
                    canonical_float(row["branch_length"]),
                )
                if replicate in timing_by_cell[timing_cell]:
                    raise ValueError(f"Duplicate timing row: {timing_cell + (replicate,)}")
                timing_by_cell[timing_cell].add(replicate)
                for field in ("simulation_seconds", "analysis_seconds", "total_seconds"):
                    value = float(row[field])
                    if not math.isfinite(value) or value < 0.0:
                        raise ValueError(f"Invalid timing value: {row}")
                timing_writer.writerow(row)
                timing_rows += 1

            if shard_index % 100 == 0 or shard_index == len(shards):
                print(
                    f"Validated shards: {shard_index}/{len(shards)}; "
                    f"truth rows: {truth_rows}; replicate rows: {replicate_rows}",
                    flush=True,
                )

        truth_handle.close()
        replicate_handle.close()
        timing_handle.close()
        truth_handle = replicate_handle = timing_handle = None
        validate_global_completeness(by_cell, expected_cells, args.expected_reps)

        timing_expected_cells = {
            (design, tree, model, str(site), f"{branch:.12g}")
            for design, tree, model, site, branch in itertools.product(
                (manifest["design"],),
                csv_values(manifest["tree_cases"]),
                csv_values(manifest["models"]),
                csv_values(manifest["site_lengths"], int),
                csv_values(manifest["branch_lengths"], float),
            )
        }
        validate_global_completeness(timing_by_cell, timing_expected_cells, args.expected_reps)
        summaries = aggregate_rows(aggregates)
        write_table(temporary["summary"], CALIBRATION_SUMMARY_FIELDS, summaries)
        heldout = []
        if manifest["design"] == "exact_null":
            heldout = heldout_calibration(training_nulls, holdout_nulls, args.alpha)
            write_table(temporary["heldout"], HELDOUT_FIELDS, heldout)

        for key in ("truth", "replicate", "summary", "timing"):
            destination = outdir / output_names[key]
            os.replace(temporary[key], destination)
            published.append(destination)
        if heldout:
            destination = outdir / output_names["heldout"]
            os.replace(temporary["heldout"], destination)
            published.append(destination)
    except Exception:
        for path in published:
            path.unlink(missing_ok=True)
        raise
    finally:
        for handle in (truth_handle, replicate_handle, timing_handle):
            if handle is not None:
                handle.close()
        for path in temporary.values():
            path.unlink(missing_ok=True)

    print(f"Shards: {len(shards)}")
    print(f"Truth rows: {truth_rows}")
    print(f"Replicate FDP rows: {replicate_rows}")
    print(f"Summary rows: {len(summaries)}")
    print(f"Timing rows: {timing_rows}")
    print(f"Held-out rows: {len(heldout)}")


if __name__ == "__main__":
    main()
