#!/usr/bin/env python3
"""Behavioral checks for the experiment discovery and execution seam."""

from __future__ import annotations

import json
import shlex
import sys
import tempfile
from pathlib import Path

from satute_analysis.experiments import ExperimentCatalog, ExperimentError
from satute_analysis.project import PROJECT_ROOT


def write_experiment(
    project: Path,
    experiment_id: str = "999_valid",
    *,
    output: str = "artifacts/local/valid",
    driver: str = "run.py",
    execution: str = "managed",
    release: str | None = None,
) -> Path:
    """Create the smallest experiment module accepted by the public interface."""

    directory = project / "experiments" / experiment_id
    directory.mkdir(parents=True)
    (directory / "README.md").write_text("# Test experiment\n", encoding="utf-8")
    (directory / "run.py").write_text(
        "if __name__ == '__main__': pass\n", encoding="utf-8"
    )
    manifest = {
        "schema_version": 1,
        "id": experiment_id,
        "question": "Does the catalog preserve its contracts?",
        "status": "test",
        "driver": driver,
        "consumers": [],
        "studies": [
            {
                "id": "contract_check",
                "question": "Does execution cross the catalog seam correctly?",
                "outcome": "diagnostic",
                "scale": "local" if execution == "managed" else "manual",
            }
        ],
        "default_study": "contract_check",
        "config": {
            "schema_version": 1,
            "output": output,
            "token": "observed",
        },
    }
    if execution == "managed":
        manifest["workflow"] = "workflow.py"
        (directory / "workflow.py").write_text(
            "def run(context, *, profile, figures, study):\n"
            "    output = context.output_dir(profile)\n"
            "    output.mkdir(parents=True, exist_ok=True)\n"
            "    token = context.experiment.configuration['token']\n"
            "    (output / 'observation.txt').write_text(\n"
            "        f'{profile}|{figures}|{study}|{token}', encoding='utf-8'\n"
            "    )\n",
            encoding="utf-8",
        )
    else:
        manifest["execution"] = execution
    if release is not None:
        manifest["release"] = release
    (directory / "experiment.json").write_text(
        json.dumps(manifest), encoding="utf-8"
    )
    return directory


def expect_discovery_error(project: Path, fragment: str) -> None:
    try:
        ExperimentCatalog.discover(project)
    except ExperimentError as error:
        if fragment not in str(error):
            raise
    else:
        raise SystemExit(f"Experiment discovery accepted invalid input: {fragment}")


def argument_value(command: list[str], option: str) -> str:
    try:
        return command[command.index(option) + 1]
    except (ValueError, IndexError) as error:
        raise SystemExit(f"Missing workflow argument {option}: {command}") from error


def check_project_catalog() -> None:
    catalog = ExperimentCatalog.discover(PROJECT_ROOT)
    numbered_directories = tuple(
        directory.name
        for directory in sorted((PROJECT_ROOT / "experiments").iterdir())
        if directory.is_dir() and directory.name[:3].isdigit()
    )
    if catalog.ids != numbered_directories:
        raise SystemExit(
            f"Automatic discovery omitted an experiment: {catalog.ids}"
        )

    for experiment in catalog:
        if experiment.manifest_path.parent != experiment.directory:
            raise SystemExit(f"Manifest is not colocated for {experiment.id}")
        if experiment.driver_path.parent != experiment.directory:
            raise SystemExit(f"Driver is not owned by {experiment.id}")
        if (
            experiment.workflow_path
            and experiment.workflow_path.parent != experiment.directory
        ):
            raise SystemExit(f"Workflow is not owned by {experiment.id}")
        if not (experiment.directory / "README.md").is_file():
            raise SystemExit(f"README is not owned by {experiment.id}")
        if (
            not experiment.question
            or not experiment.status
            or not experiment.configuration
        ):
            raise SystemExit(f"Incomplete experiment interface for {experiment.id}")
        if not experiment.studies or experiment.default_study_id not in experiment.study_ids:
            raise SystemExit(f"Incomplete study interface for {experiment.id}")

    try:
        catalog.run("005_biological_windows")
    except ExperimentError as error:
        if "requires explicit inputs" not in str(error):
            raise
    else:
        raise SystemExit("Manual biological workflow invented input defaults")

    experiment = catalog.get("002_relative_weighting")
    configuration = experiment.configuration
    configuration["schema_version"] = 999
    configuration["decision_rules"].append("invalid")
    if experiment.config_schema_version != 2:
        raise SystemExit("Experiment configuration mutated its validated schema version")
    if "invalid" in experiment.configuration["decision_rules"]:
        raise SystemExit("Experiment configuration was not returned defensively")


def check_execution_through_catalog() -> None:
    with tempfile.TemporaryDirectory(prefix="satute-experiment-catalog-") as temporary:
        project = Path(temporary)
        write_experiment(project)
        catalog = ExperimentCatalog.discover(project)
        catalog.run("999_valid", profile="smoke", figures=True, python=sys.executable)
        observation = (
            project / "artifacts/local/smoke/999_valid/schema-1/observation.txt"
        ).read_text(encoding="utf-8")
        if observation != "smoke|True|contract_check|observed":
            raise SystemExit(f"Workflow result crossed the seam incorrectly: {observation}")


def check_path_contracts() -> None:
    with tempfile.TemporaryDirectory(prefix="satute-experiment-catalog-") as temporary:
        project = Path(temporary)
        write_experiment(project, driver="../outside.py")
        expect_discovery_error(project, "escapes")

    for invalid_output in (".", "artifacts/local", "artifacts/releases/result"):
        with tempfile.TemporaryDirectory(
            prefix="satute-experiment-catalog-"
        ) as temporary:
            project = Path(temporary)
            write_experiment(project, output=invalid_output)
            expect_discovery_error(project, "dedicated directory")

    for invalid_release in (".", "artifacts/local/release", "artifacts/releases"):
        with tempfile.TemporaryDirectory(
            prefix="satute-experiment-catalog-"
        ) as temporary:
            project = Path(temporary)
            write_experiment(project, release=invalid_release)
            expect_discovery_error(project, "dedicated directory")


def check_module_completeness() -> None:
    with tempfile.TemporaryDirectory(prefix="satute-experiment-catalog-") as temporary:
        project = Path(temporary)
        directory = write_experiment(project)
        (directory / "README.md").unlink()
        expect_discovery_error(project, "README is missing")

    with tempfile.TemporaryDirectory(prefix="satute-experiment-catalog-") as temporary:
        project = Path(temporary)
        directory = write_experiment(project)
        (directory / "experiment.json").unlink()
        expect_discovery_error(project, "missing experiment.json")


def check_unique_outputs() -> None:
    with tempfile.TemporaryDirectory(prefix="satute-experiment-catalog-") as temporary:
        project = Path(temporary)
        write_experiment(
            project, "998_first", output="artifacts/local/shared"
        )
        write_experiment(
            project, "999_second", output="artifacts/local/shared"
        )
        expect_discovery_error(project, "Duplicate experiment output directories")


def check_relative_weighting_workflow_contract() -> None:
    catalog = ExperimentCatalog.discover(PROJECT_ROOT)
    experiment = catalog.get("002_relative_weighting")
    design = experiment.configuration["publication_curve_design"]
    with tempfile.TemporaryDirectory(prefix="satute-workflow-command-") as temporary:
        temporary_path = Path(temporary)
        command_log = temporary_path / "command.txt"
        fake_python = temporary_path / "python"
        fake_python.write_text(
            "#!/bin/sh\n"
            f"printf '%s\\n' \"$@\" > {shlex.quote(str(command_log))}\n",
            encoding="utf-8",
        )
        fake_python.chmod(0o755)

        expected_power_studies = {
            "gtr_power",
            "protein_power",
            "misspecification_power",
        }
        if set(design["studies"]) != expected_power_studies:
            raise SystemExit("A declared relative-weighting power study is missing")
        expected_catalog_studies = expected_power_studies | {
            "fdr_exact_null",
            "fdr_mixed_null",
        }
        if set(experiment.study_ids) != expected_catalog_studies:
            raise SystemExit("A power or calibration study disappeared from the catalog")
        for study_id, study_design in design["studies"].items():
            catalog.run(
                experiment.id,
                study=study_id,
                profile="full",
                python=fake_python,
            )
            command = command_log.read_text(encoding="utf-8").splitlines()
            expected = {
                "--reps": str(design["replicates_per_cell"]),
                "--site-lengths": ",".join(map(str, design["site_lengths"])),
                "--branch-lengths": ",".join(map(str, design["branch_lengths"])),
                "--tree-cases": ",".join(design["tree_cases"]),
                "--model-pairs": ",".join(study_design["model_pairs"]),
                "--scenario-set": study_design["scenario_set"],
                "--simulator": design["simulator"],
            }
            for option, value in expected.items():
                if argument_value(command, option) != value:
                    raise SystemExit(f"{study_id} does not honor manifest {option}")
            if "--paper-grid" in command:
                raise SystemExit("Full workflow silently replaced the manifest grid")

        calibration = experiment.configuration["fdr_calibration_design"]
        for study_id, driver_design, expected_reps in (
            (
                "fdr_exact_null",
                "exact_null",
                calibration["exact_null_replicates_per_cell"],
            ),
            (
                "fdr_mixed_null",
                "mixed_null",
                calibration["mixed_null_replicates_per_cell"],
            ),
        ):
            catalog.run(
                experiment.id,
                study=study_id,
                profile="full",
                python=fake_python,
            )
            command = command_log.read_text(encoding="utf-8").splitlines()
            if not command[0].endswith("fdr_calibration_run.py"):
                raise SystemExit(f"{study_id} lost its calibration driver")
            if argument_value(command, "--design") != driver_design:
                raise SystemExit(f"{study_id} selected the wrong calibration design")
            if argument_value(command, "--reps") != str(expected_reps):
                raise SystemExit(f"{study_id} lost its declared replicate count")

        catalog.run(
            experiment.id,
            study="gtr_power",
            profile="smoke",
            python=fake_python,
        )
        smoke_command = command_log.read_text(encoding="utf-8").splitlines()
        if "--resume" in smoke_command:
            raise SystemExit("Smoke workflow can reuse stale output")


def main() -> None:
    check_project_catalog()
    check_execution_through_catalog()
    check_path_contracts()
    check_module_completeness()
    check_unique_outputs()
    check_relative_weighting_workflow_contract()
    print("Experiment catalog contracts: passed")


if __name__ == "__main__":
    main()
