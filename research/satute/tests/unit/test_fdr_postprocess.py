#!/usr/bin/env python3
"""End-to-end gates for the streaming FDR calibration postprocessor."""

from __future__ import annotations

import csv
import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from satute_analysis.benchmark_contracts import SCHEMA_VERSION
from satute_analysis.fdr_calibration import (
    CALIBRATION_SUMMARY_FIELDS,
    REPLICATE_FDP_FIELDS,
    TRUTH_BRANCH_FIELDS,
)


PROJECT = Path(__file__).resolve().parents[2]
POSTPROCESS = PROJECT / "experiments/002_relative_weighting/fdr_calibration_postprocess.py"
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


def load_postprocess_module():
    specification = importlib.util.spec_from_file_location(
        "satute_fdr_calibration_postprocess", POSTPROCESS
    )
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def check_discrete_tie_calibration():
    module = load_postprocess_module()
    cell = ("exact_null", "five_external", "JC", "JC", "100", "4", "dominant", "taxon_bonferroni")

    # A mass of adjusted p-values at one must never turn a nominal 5% rule
    # into 100% rejection simply because the selected order statistic is tied.
    all_tied = module.heldout_calibration({cell: [1.0] * 20}, {cell: [1.0] * 20}, 0.05)[0]
    if not all_tied["empirical_threshold"] < 1.0:
        raise SystemExit("Discrete all-tie calibration did not move below the minimum")
    if all_tied["holdout_rejections"] != 0:
        raise SystemExit("Discrete all-tie calibration is anti-conservative")

    one_attainable = module.heldout_calibration(
        {cell: [0.01] + [1.0] * 19},
        {cell: [0.01, 0.02] + [1.0] * 18},
        0.05,
    )[0]
    if one_attainable["empirical_threshold"] != 0.01:
        raise SystemExit("Calibration failed to retain the largest valid threshold")
    if one_attainable["holdout_rejections"] != 1:
        raise SystemExit("Held-out rejection count does not match the calibrated rule")


def write_table(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def fixture(root):
    shard = root / "shards/shard-0"
    manifest = {
        "schema_version": "2",
        "design": "exact_null",
        "replicates": "2",
        "training_replicates": "1",
        "shard_count": "1",
        "models": "JC",
        "resolved_model_pairs": "JC=>JC",
        "site_lengths": "100",
        "branch_lengths": "4",
        "tree_cases": "five_external",
    }
    write_table(
        root / "run_manifest.tsv",
        ("field", "value"),
        ({"field": key, "value": value} for key, value in manifest.items()),
    )
    replicate_rows = []
    truth_rows = []
    timing_rows = []
    for replicate in (1, 2):
        sample = "training" if replicate == 1 else "holdout"
        seed = 100 + replicate
        for formula in ("dominant", "eigenvalue_weighted"):
            for rule in ("unadjusted", "taxon_bonferroni", "by_fdr"):
                base = {
                    "schema_version": SCHEMA_VERSION,
                    "design": "exact_null",
                    "sample": sample,
                    "tree_case": "five_external",
                    "simulation_model": "JC",
                    "evaluation_model": "JC",
                    "nsites": 100,
                    "branch_length": 4.0,
                    "replicate": replicate,
                    "seed": seed,
                    "formula": formula,
                    "decision_rule": rule,
                }
                replicate_rows.append(
                    {
                        **base,
                        "family_size": 2,
                        "null_branches": 1,
                        "alternative_branches": 1,
                        "discoveries": 1,
                        "false_discoveries": 0,
                        "true_discoveries": 1,
                        "fdp": 0.0,
                        "power": 1.0,
                    }
                )
                for branch_id, split, truth, p_value, decision in (
                    ("1", "B", "null", 0.8, "saturated"),
                    ("2", "A", "alternative", 0.01, "informative"),
                ):
                    discovery = int(decision == "informative")
                    truth_rows.append(
                        {
                            **base,
                            "branch_id": branch_id,
                            "split": split,
                            "truth": truth,
                            "p_value": p_value,
                            "decision": decision,
                            "discovery": discovery,
                            "false_discovery": int(discovery and truth == "null"),
                            "true_discovery": int(discovery and truth == "alternative"),
                        }
                    )
        timing_rows.append(
            {
                "schema_version": SCHEMA_VERSION,
                "design": "exact_null",
                "tree_case": "five_external",
                "simulation_model": "JC",
                "nsites": 100,
                "branch_length": 4.0,
                "replicate": replicate,
                "seed": seed,
                "truth": "null",
                "simulation_seconds": 0.01,
                "analysis_seconds": 0.02,
                "total_seconds": 0.03,
            }
        )
    write_table(shard / "replicate_fdp.tsv", REPLICATE_FDP_FIELDS, replicate_rows)
    write_table(shard / "truth_branch_audit.tsv", TRUTH_BRANCH_FIELDS, truth_rows)
    write_table(shard / "calibration_summary.tsv", CALIBRATION_SUMMARY_FIELDS, [])
    write_table(shard / "timing.tsv", TIMING_FIELDS, timing_rows)


def run(root, outdir, expect_success):
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(PROJECT / "src") + os.pathsep + environment.get(
        "PYTHONPATH", ""
    )
    result = subprocess.run(
        [
            sys.executable,
            str(POSTPROCESS),
            "--root",
            str(root),
            "--outdir",
            str(outdir),
            "--expected-shards",
            "1",
            "--expected-reps",
            "2",
            "--training-reps",
            "1",
        ],
        cwd=PROJECT,
        env=environment,
        capture_output=True,
        text=True,
    )
    if expect_success and result.returncode:
        raise SystemExit(f"Streaming postprocessor failed:\n{result.stdout}\n{result.stderr}")
    if not expect_success and not result.returncode:
        raise SystemExit("Corrupt truth audit unexpectedly passed postprocessing")
    return result


def main():
    check_discrete_tie_calibration()
    with tempfile.TemporaryDirectory(prefix="satute-fdr-postprocess-") as temporary:
        temporary = Path(temporary)
        root = temporary / "valid"
        fixture(root)
        outdir = temporary / "valid-output"
        run(root, outdir, expect_success=True)
        expected_outputs = {
            "truth_branch_audit.tsv",
            "replicate_fdp.tsv",
            "calibration_summary.tsv",
            "timing.tsv",
            "heldout_null_calibration.tsv",
        }
        if {path.name for path in outdir.iterdir()} != expected_outputs:
            raise SystemExit(f"Unexpected streaming outputs: {sorted(outdir.iterdir())}")
        with (outdir / "replicate_fdp.tsv").open(encoding="utf-8") as handle:
            if sum(1 for _line in handle) - 1 != 12:
                raise SystemExit("Merged replicate row count is incorrect")

        corrupt = temporary / "corrupt"
        shutil.copytree(root, corrupt)
        truth_path = corrupt / "shards/shard-0/truth_branch_audit.tsv"
        with truth_path.open("r", encoding="utf-8") as handle:
            lines = handle.readlines()
        with truth_path.open("w", encoding="utf-8") as handle:
            handle.writelines(lines[:-1])
        corrupt_output = temporary / "corrupt-output"
        result = run(corrupt, corrupt_output, expect_success=False)
        if "Truth-family" not in result.stderr:
            raise SystemExit(f"Corruption failed for the wrong reason: {result.stderr}")
        if corrupt_output.exists() and any(corrupt_output.iterdir()):
            raise SystemExit("Failed postprocessing left accepted-looking output files")
    print("Streaming FDR postprocessing checks passed")


if __name__ == "__main__":
    main()
