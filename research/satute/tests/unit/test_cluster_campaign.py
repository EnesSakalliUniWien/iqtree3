#!/usr/bin/env python3
"""Contract checks for the versioned LiSC recalculation campaign."""

from __future__ import annotations

import subprocess

from satute_analysis.cluster_campaign import (
    DEFAULT_CAMPAIGN,
    branch_lengths,
    load_campaign,
    task_count,
)
from satute_analysis.experiments import ExperimentCatalog
from satute_analysis.project import PROJECT_ROOT


def main() -> None:
    campaign = load_campaign()
    runs = {run["id"]: run for run in campaign["runs"]}
    expected = {
        "jc_power",
        "gtr_power",
        "protein_power",
        "misspecification_power",
        "fdr_exact_null",
        "fdr_mixed_null",
        "corrected_invariant_mixture",
        "saturation_diagnostics",
        "biological_windows",
    }
    if set(runs) != expected:
        raise SystemExit(f"Campaign coverage changed: {sorted(runs)}")
    precision = campaign["precision_contract"]
    if precision["main_power_replicates_per_cell"] != 10000:
        raise SystemExit("Campaign-level power precision disagrees with its runs")
    if precision["exact_null_holdout_replicates"] != 10000:
        raise SystemExit("Campaign-level null precision target changed")

    if runs["jc_power"]["baseline_branch_lengths"][-1] != 10:
        raise SystemExit("JC original-paper baseline was not preserved")
    for run_id in ("jc_power", "gtr_power", "protein_power", "misspecification_power"):
        run = runs[run_id]
        if run["replicates_per_cell"] != 10000:
            raise SystemExit(f"{run_id} lost the precision target")
        if max(branch_lengths(run)) < 50:
            raise SystemExit(f"{run_id} lost the long-branch extension")
        if task_count(run) is None or task_count(run) <= 0:
            raise SystemExit(f"{run_id} has no computable task count")

    exact = runs["fdr_exact_null"]
    if exact["replicates_per_cell"] != 20000:
        raise SystemExit("Exact-null calibration lost its precision target")
    if exact["training_replicates"] * 2 != exact["replicates_per_cell"]:
        raise SystemExit("Exact-null training/holdout split is not balanced")
    if runs["fdr_mixed_null"]["replicates_per_cell"] != 10000:
        raise SystemExit("Mixed-null calibration lost its precision target")
    if runs["corrected_invariant_mixture"]["replicates_per_cell"] != 40000:
        raise SystemExit("Invariant-mixture precision was silently reduced")
    catalog = ExperimentCatalog.discover(PROJECT_ROOT)
    invariant = catalog.get("003_invariant_mixture").configuration
    if runs["corrected_invariant_mixture"]["replicates_per_cell"] != (
        invariant["null_replicates"] + invariant["alternative_replicates"]
    ):
        raise SystemExit("Campaign and invariant-mixture manifest disagree")
    if runs["corrected_invariant_mixture"]["native_iqtree_equivalent"]:
        raise SystemExit("Independent correction was mislabeled as native IQ-TREE")
    saturation = catalog.get("004_saturation_diagnostics").configuration
    if runs["saturation_diagnostics"]["replicates_per_cell"] != saturation[
        "alternative_replicates"
    ]:
        raise SystemExit("Campaign and saturation manifest disagree on alternatives")
    if runs["saturation_diagnostics"]["null_replicates"] != saturation[
        "null_replicates"
    ]:
        raise SystemExit("Campaign and saturation manifest disagree on nulls")
    if runs["biological_windows"]["execution"] != "manual":
        raise SystemExit("Biological input requirements were hidden")

    result = subprocess.run(
        [
            str(PROJECT_ROOT / "satute"),
            "cluster",
            "plan",
            "--campaign",
            str(DEFAULT_CAMPAIGN),
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    for run_id in expected:
        if f"{run_id}:" not in result.stdout:
            raise SystemExit(f"Cluster plan omitted {run_id}")

    print("Expanded LiSC campaign contract: passed")


if __name__ == "__main__":
    main()
