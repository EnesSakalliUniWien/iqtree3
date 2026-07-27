"""Plan and submit versioned SatuTe recalculation campaigns on LiSC."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .project import PROJECT_ROOT


DEFAULT_CAMPAIGN = PROJECT_ROOT / "campaigns/2026_expanded_recalculation.json"
LISC_TOOLS = PROJECT_ROOT / "tools/cluster/lisc"
REPOSITORY_ROOT = PROJECT_ROOT.parents[1]


class CampaignError(RuntimeError):
    """Raised when a campaign cannot be planned or submitted safely."""


def load_campaign(path: Path = DEFAULT_CAMPAIGN) -> dict[str, Any]:
    campaign = json.loads(path.read_text(encoding="utf-8"))
    if campaign.get("schema_version") != 1:
        raise CampaignError("Unsupported campaign schema")
    runs = campaign.get("runs")
    if not isinstance(runs, list) or not runs:
        raise CampaignError("Campaign must declare at least one run")
    ids = [run.get("id") for run in runs]
    if any(not run_id for run_id in ids) or len(ids) != len(set(ids)):
        raise CampaignError("Campaign run IDs must be present and unique")
    return campaign


def branch_lengths(run: dict[str, Any]) -> list[float]:
    values = [
        *run.get("baseline_branch_lengths", []),
        *run.get("extension_branch_lengths", []),
    ]
    if len(values) != len(set(map(float, values))):
        raise CampaignError(f"{run['id']} repeats a branch length")
    return values


def task_count(run: dict[str, Any]) -> int | None:
    if "replicates_per_cell" not in run:
        return None
    dimensions = (
        len(run.get("tree_cases", [])),
        len(run.get("site_lengths", [])),
        len(branch_lengths(run)),
        len(run.get("model_pairs", run.get("models", []))),
        int(run["replicates_per_cell"]),
    )
    total = 1
    for size in dimensions:
        total *= size
    return total


def selected_runs(
    campaign: dict[str, Any], requested: list[str] | None = None
) -> list[dict[str, Any]]:
    runs = campaign["runs"]
    if not requested:
        return runs
    by_id = {run["id"]: run for run in runs}
    unknown = sorted(set(requested) - by_id.keys())
    if unknown:
        raise CampaignError(f"Unknown campaign runs: {unknown}")
    return [by_id[run_id] for run_id in requested]


def print_plan(campaign: dict[str, Any], requested: list[str] | None = None) -> None:
    print(f"Campaign: {campaign['id']}")
    print(campaign["description"])
    precision = campaign.get("precision_contract", {})
    if precision:
        print(
            "Precision: "
            f"{int(precision['main_power_replicates_per_cell']):,} power "
            "replicates/cell; worst-case 95% half-width "
            f"{100 * float(precision['main_power_worst_case_95_percent_half_width']):.2f} pp"
        )
    resources = campaign.get("resource_estimate", {})
    if resources:
        print(
            "Estimated cluster scale: "
            f"{int(resources['cluster_cpu_hours']):,} CPU-hours; "
            f"{int(resources['cluster_scratch_gib']):,} GiB scratch"
        )
    print()
    for run in selected_runs(campaign, requested):
        baseline = run.get("baseline_branch_lengths", [])
        extension = run.get("extension_branch_lengths", [])
        tasks = task_count(run)
        print(f"{run['id']}: {run['experiment']}/{run['study']}")
        print(f"  execution: {run['execution']}")
        if "replicates_per_cell" in run:
            print(f"  replicates/cell: {int(run['replicates_per_cell']):,}")
        if tasks is not None:
            print(f"  replicate tasks: {tasks:,}")
        if baseline:
            print("  baseline branches: " + ",".join(map(str, baseline)))
        if extension:
            print("  extension branches: " + ",".join(map(str, extension)))
        if note := run.get("replicate_note"):
            print(f"  note: {note}")
        if command := run.get("command"):
            print("  command: " + " ".join(command))
        if command := run.get("postprocess_command"):
            print("  postprocess: " + " ".join(command))
        print()


def _git_output(*arguments: str) -> str:
    result = subprocess.run(
        ["git", "-C", str(REPOSITORY_ROOT), *arguments],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def _submission_environment(
    campaign: dict[str, Any], run: dict[str, Any], source_id: str
) -> dict[str, str]:
    cluster = campaign["cluster"]
    branches = branch_lengths(run)
    environment = dict(os.environ)
    environment.update(
        {
            "ROOT_DIR": str(REPOSITORY_ROOT),
            "PROJECT_DIR": str(PROJECT_ROOT),
            "WORK_DIR": cluster["work_dir"],
            "SOURCE_ID": source_id,
            "REPS": str(run["replicates_per_cell"]),
            "SHARD_COUNT": str(cluster["shard_count"]),
            "TREE_CASES": ",".join(run["tree_cases"]),
            "SITE_LENGTHS": ",".join(map(str, run["site_lengths"])),
            "BRANCH_LENGTHS": ",".join(map(str, branches)),
            "RUN_NAME": (
                f"{campaign['id']}-{run['id']}-reps"
                f"{run['replicates_per_cell']}-{source_id}"
            ),
        }
    )
    if run["submitter"] == "fdr":
        environment["MODELS"] = ",".join(run["models"])
        environment["TRAINING_REPS"] = str(run["training_replicates"])
    return environment


def _submission_command(run: dict[str, Any]) -> list[str]:
    slurm = run["slurm"]
    options = [
        f"--array=0-999%{slurm['array_limit']}",
        f"--time={slurm['time']}",
        f"--mem={slurm['memory']}",
    ]
    if run["submitter"] == "publication":
        script = LISC_TOOLS / "submit_publication_suite.sh"
    elif run["submitter"] == "fdr":
        script = LISC_TOOLS / "submit_fdr_calibration_suite.sh"
    else:
        raise CampaignError(f"Unknown submitter for {run['id']}")
    return ["bash", str(script), run["design"], *options]


def submit_campaign(
    campaign: dict[str, Any], requested: list[str] | None = None
) -> Path:
    runs = [
        run
        for run in selected_runs(campaign, requested)
        if run["execution"] == "cluster"
    ]
    if not runs:
        raise CampaignError("Selection contains no cluster runs")
    if shutil.which("sbatch") is None:
        raise CampaignError(
            "sbatch is unavailable; run submission from a LiSC login node"
        )
    if campaign["source_policy"]["require_clean_tracked_tree"]:
        dirty = _git_output("status", "--short", "--untracked-files=no")
        if dirty:
            raise CampaignError(
                "Tracked files are modified; commit the exact scientific source "
                "before cluster submission"
            )
    source_id = _git_output("rev-parse", "--short=12", "HEAD")
    iqtree = (
        Path(campaign["cluster"]["work_dir"])
        / f"iqtree3-build-{source_id}"
        / "iqtree3"
    )
    if (
        campaign["source_policy"]["require_source_matched_iqtree"]
        and not os.access(iqtree, os.X_OK)
    ):
        raise CampaignError(f"Missing source-matched IQ-TREE binary: {iqtree}")

    records = []
    for run in runs:
        environment = _submission_environment(campaign, run, source_id)
        command = _submission_command(run)
        result = subprocess.run(
            command,
            cwd=REPOSITORY_ROOT,
            env=environment,
            check=True,
            capture_output=True,
            text=True,
        )
        job_id = result.stdout.strip()
        if not job_id.isdigit():
            raise CampaignError(
                f"{run['id']} returned an invalid Slurm job ID: {job_id!r}"
            )
        records.append(
            {
                "run_id": run["id"],
                "job_id": job_id,
                "run_name": environment["RUN_NAME"],
                "command": command,
            }
        )
        print(f"{run['id']}: job {job_id}")

    submitted = datetime.now(timezone.utc)
    output = (
        PROJECT_ROOT
        / "artifacts/cluster/lisc/submissions"
        / f"{campaign['id']}-{submitted:%Y%m%dT%H%M%SZ}"
    )
    output.mkdir(parents=True, exist_ok=False)
    manifest = {
        "schema_version": 1,
        "campaign": campaign["id"],
        "submitted_utc": submitted.isoformat(),
        "source_id": source_id,
        "iqtree": str(iqtree),
        "runs": records,
    }
    path = output / "submission_manifest.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(path)
    return path


__all__ = [
    "CampaignError",
    "DEFAULT_CAMPAIGN",
    "branch_lengths",
    "load_campaign",
    "print_plan",
    "selected_runs",
    "submit_campaign",
    "task_count",
]
