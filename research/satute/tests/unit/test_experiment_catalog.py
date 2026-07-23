#!/usr/bin/env python3
"""Behavioral checks for the experiment discovery and execution seam."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

from satute_analysis.experiments import (
    ExperimentCatalog,
    ExperimentError,
    WorkflowContext,
)
from satute_analysis.project import PROJECT_ROOT


EXPECTED_IDS = (
    "001_paper_reconstruction",
    "002_relative_weighting",
    "003_invariant_mixture",
    "004_saturation_diagnostics",
    "005_biological_windows",
)


def check_project_catalog() -> None:
    catalog = ExperimentCatalog.discover(PROJECT_ROOT)
    if catalog.ids != EXPECTED_IDS:
        raise SystemExit(f"Unexpected experiment catalog: {catalog.ids}")

    for experiment in catalog:
        if experiment.manifest_path.parent != experiment.directory:
            raise SystemExit(f"Manifest is not colocated for {experiment.id}")
        if experiment.driver_path.parent != experiment.directory:
            raise SystemExit(f"Driver is not owned by {experiment.id}")
        if experiment.workflow_path.parent != experiment.directory:
            raise SystemExit(f"Workflow is not owned by {experiment.id}")
        if not experiment.question or not experiment.status or not experiment.configuration:
            raise SystemExit(f"Incomplete experiment interface for {experiment.id}")
        context = WorkflowContext(PROJECT_ROOT, experiment, Path("python3"))
        expected_output = (
            PROJECT_ROOT
            / "artifacts/local/smoke"
            / experiment.id
            / f"schema-{experiment.config_schema_version}"
        )
        if context.output_dir("smoke") != expected_output:
            raise SystemExit(f"Smoke output is not schema-isolated for {experiment.id}")

    try:
        catalog.run("005_biological_windows")
    except ExperimentError as error:
        if "explicit alignment inputs" not in str(error):
            raise
    else:
        raise SystemExit("Biological workflow invented defaults instead of requiring inputs")

    experiment = catalog.get("002_relative_weighting")
    configuration = experiment.configuration
    configuration["schema_version"] = 999
    configuration["decision_rules"].append("invalid")
    if experiment.config_schema_version != 2:
        raise SystemExit("Experiment configuration mutated its validated schema version")
    if "invalid" in experiment.configuration["decision_rules"]:
        raise SystemExit("Experiment configuration was not returned defensively")


def check_path_containment() -> None:
    with tempfile.TemporaryDirectory(prefix="satute-experiment-catalog-") as temporary:
        project = Path(temporary)
        experiment_dir = project / "experiments/999_invalid"
        experiment_dir.mkdir(parents=True)
        (project / "outside.py").write_text("def run(*args, **kwargs): pass\n")
        (experiment_dir / "workflow.py").write_text("def run(*args, **kwargs): pass\n")
        manifest = {
            "schema_version": 1,
            "id": "999_invalid",
            "question": "Can an owned path escape its experiment module?",
            "status": "test",
            "driver": "../outside.py",
            "workflow": "workflow.py",
            "consumers": [],
            "config": {"output": "artifacts/local/invalid"},
        }
        (experiment_dir / "experiment.json").write_text(json.dumps(manifest))
        try:
            ExperimentCatalog.discover(project)
        except ExperimentError as error:
            if "escapes its experiment module" not in str(error):
                raise
        else:
            raise SystemExit("Experiment path containment was not enforced")


def check_project_path_containment() -> None:
    with tempfile.TemporaryDirectory(prefix="satute-experiment-catalog-") as temporary:
        project = Path(temporary)
        experiment_dir = project / "experiments/999_invalid"
        experiment_dir.mkdir(parents=True)
        (experiment_dir / "run.py").write_text("if __name__ == '__main__': pass\n")
        (experiment_dir / "workflow.py").write_text("def run(*args, **kwargs): pass\n")
        manifest = {
            "schema_version": 1,
            "id": "999_invalid",
            "question": "Can an output path escape the research project?",
            "status": "test",
            "driver": "run.py",
            "workflow": "workflow.py",
            "consumers": [],
            "config": {
                "schema_version": 1,
                "output": "../outside",
            },
        }
        (experiment_dir / "experiment.json").write_text(json.dumps(manifest))
        try:
            ExperimentCatalog.discover(project)
        except ExperimentError as error:
            if "escapes the SatuTe project" not in str(error):
                raise
        else:
            raise SystemExit("Experiment output path containment was not enforced")


def check_unique_outputs() -> None:
    with tempfile.TemporaryDirectory(prefix="satute-experiment-catalog-") as temporary:
        project = Path(temporary)
        for experiment_id in ("998_first", "999_second"):
            experiment_dir = project / f"experiments/{experiment_id}"
            experiment_dir.mkdir(parents=True)
            (experiment_dir / "run.py").write_text("if __name__ == '__main__': pass\n")
            (experiment_dir / "workflow.py").write_text(
                "def run(*args, **kwargs): pass\n"
            )
            manifest = {
                "schema_version": 1,
                "id": experiment_id,
                "question": "Can two experiments claim one output directory?",
                "status": "test",
                "driver": "run.py",
                "workflow": "workflow.py",
                "consumers": [],
                "config": {
                    "schema_version": 1,
                    "output": "artifacts/local/shared",
                },
            }
            (experiment_dir / "experiment.json").write_text(json.dumps(manifest))
        try:
            ExperimentCatalog.discover(project)
        except ExperimentError as error:
            if "Duplicate experiment output directories" not in str(error):
                raise
        else:
            raise SystemExit("Duplicate experiment output directories were accepted")


def main() -> None:
    check_project_catalog()
    check_path_containment()
    check_project_path_containment()
    check_unique_outputs()
    print("Experiment catalog contracts: passed")


if __name__ == "__main__":
    main()
