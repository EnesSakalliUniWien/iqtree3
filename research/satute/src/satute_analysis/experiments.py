"""Discovery and execution interface for the SatuTe experiment modules."""

from __future__ import annotations

import importlib.util
import json
import os
import re
import shlex
import subprocess
import sys
from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterator, Mapping

from .project import PROJECT_ROOT


MANIFEST_NAME = "experiment.json"
EXPERIMENT_ID = re.compile(r"^\d{3}_[a-z0-9]+(?:_[a-z0-9]+)*$")
REQUIRED_FIELDS = {
    "schema_version",
    "id",
    "question",
    "status",
    "driver",
    "workflow",
    "consumers",
    "config",
}


class ExperimentError(RuntimeError):
    """Raised when an experiment definition or workflow is invalid."""


@dataclass(frozen=True)
class Experiment:
    """Validated, immutable description of one experiment module."""

    directory: Path
    id: str
    question: str
    status: str
    consumers: tuple[str, ...]
    driver_path: Path
    workflow_path: Path
    output_path: Path
    config_schema_version: int
    release_path: Path | None
    _configuration: dict[str, Any] = field(repr=False)

    @property
    def manifest_path(self) -> Path:
        return self.directory / MANIFEST_NAME

    @property
    def configuration(self) -> dict[str, Any]:
        """Return a defensive copy suitable for provenance serialization."""

        return deepcopy(self._configuration)

    def file(self, relative_path: str | Path) -> Path:
        """Resolve a file owned by this experiment without escaping its module."""

        if Path(relative_path).is_absolute():
            raise ExperimentError(
                f"{self.id} path must be relative to its experiment module: {relative_path}"
            )
        path = (self.directory / relative_path).resolve()
        try:
            path.relative_to(self.directory.resolve())
        except ValueError as error:
            raise ExperimentError(
                f"{self.id} path escapes its experiment module: {relative_path}"
            ) from error
        return path


@dataclass(frozen=True)
class WorkflowContext:
    """Capabilities supplied to experiment workflow adapters."""

    project_root: Path
    experiment: Experiment
    python: Path

    @property
    def repository_root(self) -> Path:
        return self.project_root.parents[1]

    @property
    def iqtree(self) -> Path:
        return self.repository_root / "build/iqtree3"

    def output_dir(self, profile: str) -> Path:
        if profile == "smoke":
            return (
                self.project_root
                / "artifacts/local/smoke"
                / self.experiment.id
                / f"schema-{self.experiment.config_schema_version}"
            )
        return self.experiment.output_path

    def run_python(
        self, script: str | Path, *arguments: object
    ) -> subprocess.CompletedProcess[bytes]:
        """Run an owned Python driver with the project import path configured."""

        command = [
            str(self.python),
            str(self.experiment.file(script)),
            *(str(argument) for argument in arguments),
        ]
        environment = dict(os.environ)
        environment["PYTHONPATH"] = str(self.project_root / "src") + os.pathsep + environment.get(
            "PYTHONPATH", ""
        )
        print("+", shlex.join(command))
        return subprocess.run(
            command,
            cwd=self.project_root,
            env=environment,
            check=True,
        )


class ExperimentCatalog:
    """Validated collection and execution seam for all experiment modules."""

    def __init__(self, project_root: Path, experiments: tuple[Experiment, ...]):
        self.project_root = project_root
        self._experiments = experiments
        self._by_id = {experiment.id: experiment for experiment in experiments}

    @classmethod
    def discover(cls, project_root: Path = PROJECT_ROOT) -> "ExperimentCatalog":
        project_root = Path(project_root).resolve()
        experiment_root = project_root / "experiments"
        experiments = tuple(
            cls._read_manifest(path)
            for path in sorted(experiment_root.glob(f"*/{MANIFEST_NAME}"))
        )
        if not experiments:
            raise ExperimentError(f"No experiment manifests found under {experiment_root}")
        ids = [experiment.id for experiment in experiments]
        if len(ids) != len(set(ids)):
            raise ExperimentError(f"Duplicate experiment ids: {ids}")
        outputs = [experiment.output_path for experiment in experiments]
        if len(outputs) != len(set(outputs)):
            raise ExperimentError(f"Duplicate experiment output directories: {outputs}")
        return cls(project_root, experiments)

    @classmethod
    def _read_manifest(cls, path: Path) -> Experiment:
        try:
            manifest = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise ExperimentError(f"Cannot read experiment manifest {path}: {error}") from error

        if not isinstance(manifest, dict):
            raise ExperimentError(f"{path} must contain a JSON object")
        missing = sorted(REQUIRED_FIELDS - manifest.keys())
        if missing:
            raise ExperimentError(f"{path} is missing required fields: {missing}")

        experiment_id = cls._required_text(manifest, "id", path)
        if not EXPERIMENT_ID.fullmatch(experiment_id):
            raise ExperimentError(f"Invalid experiment id: {experiment_id}")
        if path.parent.name != experiment_id:
            raise ExperimentError(
                f"Experiment id {experiment_id} does not match directory {path.parent.name}"
            )
        if type(manifest["schema_version"]) is not int or manifest["schema_version"] != 1:
            raise ExperimentError(
                f"Unsupported manifest schema for {experiment_id}: "
                f"{manifest['schema_version']}"
            )

        question = cls._required_text(manifest, "question", path)
        status = cls._required_text(manifest, "status", path)
        driver = cls._required_text(manifest, "driver", path)
        workflow = cls._required_text(manifest, "workflow", path)

        if not isinstance(manifest["consumers"], list):
            raise ExperimentError(f"{experiment_id} consumers must be a list")
        consumers = tuple(manifest["consumers"])
        if any(not isinstance(consumer, str) or not consumer.strip() for consumer in consumers):
            raise ExperimentError(f"{experiment_id} consumers must contain non-empty strings")
        if len(consumers) != len(set(consumers)):
            raise ExperimentError(f"{experiment_id} consumers must be unique")

        configuration = manifest["config"]
        if not isinstance(configuration, dict):
            raise ExperimentError(f"{experiment_id} config must be an object")
        config_schema_version = configuration.get("schema_version", 1)
        if type(config_schema_version) is not int or config_schema_version < 1:
            raise ExperimentError(
                f"{experiment_id} config schema_version must be a positive integer"
            )
        output = cls._required_text(configuration, "output", path)
        project_root = path.parent.parent.parent.resolve()
        output_path = cls._project_path(
            project_root, output, f"{experiment_id} output directory"
        )

        release_value = manifest.get("release")
        if release_value is not None and (
            not isinstance(release_value, str) or not release_value.strip()
        ):
            raise ExperimentError(f"{experiment_id} release must be a non-empty string")
        release_path = (
            cls._project_path(project_root, release_value, f"{experiment_id} release")
            if release_value
            else None
        )
        directory = path.parent.resolve()

        experiment = Experiment(
            directory=directory,
            id=experiment_id,
            question=question,
            status=status,
            consumers=consumers,
            driver_path=cls._owned_path(directory, driver, experiment_id, "driver"),
            workflow_path=cls._owned_path(directory, workflow, experiment_id, "workflow"),
            output_path=output_path,
            config_schema_version=config_schema_version,
            release_path=release_path,
            _configuration=deepcopy(configuration),
        )
        for role, owned_path in (
            ("driver", experiment.driver_path),
            ("workflow", experiment.workflow_path),
        ):
            if not owned_path.is_file():
                raise ExperimentError(f"{experiment.id} {role} is missing: {owned_path}")
        return experiment

    @staticmethod
    def _required_text(
        values: Mapping[str, Any], field_name: str, source: Path
    ) -> str:
        value = values.get(field_name)
        if not isinstance(value, str) or not value.strip():
            raise ExperimentError(f"{source} field {field_name!r} must be a non-empty string")
        return value

    @staticmethod
    def _project_path(project_root: Path, relative_path: str, label: str) -> Path:
        path = Path(relative_path)
        if path.is_absolute():
            raise ExperimentError(f"{label} must be relative to {project_root}: {relative_path}")
        resolved = (project_root / path).resolve()
        try:
            resolved.relative_to(project_root)
        except ValueError as error:
            raise ExperimentError(
                f"{label} escapes the SatuTe project: {relative_path}"
            ) from error
        return resolved

    @staticmethod
    def _owned_path(
        directory: Path, relative_path: str, experiment_id: str, role: str
    ) -> Path:
        path = Path(relative_path)
        if path.is_absolute():
            raise ExperimentError(
                f"{experiment_id} {role} path must be relative: {relative_path}"
            )
        resolved = (directory / path).resolve()
        try:
            resolved.relative_to(directory)
        except ValueError as error:
            raise ExperimentError(
                f"{experiment_id} path escapes its experiment module: {relative_path}"
            ) from error
        return resolved

    def __iter__(self) -> Iterator[Experiment]:
        return iter(self._experiments)

    @property
    def ids(self) -> tuple[str, ...]:
        return tuple(experiment.id for experiment in self._experiments)

    def get(self, experiment_id: str) -> Experiment:
        try:
            return self._by_id[experiment_id]
        except KeyError as error:
            raise ExperimentError(
                f"Unknown experiment {experiment_id}; choose from {list(self.ids)}"
            ) from error

    def run(
        self,
        experiment_id: str,
        *,
        profile: str = "smoke",
        figures: bool = False,
        python: Path | str = sys.executable,
    ) -> None:
        if profile not in {"smoke", "full"}:
            raise ExperimentError(f"Unknown profile {profile!r}; choose smoke or full")
        experiment = self.get(experiment_id)
        workflow = self._load_workflow(experiment)
        context = WorkflowContext(self.project_root, experiment, Path(python))
        workflow(context, profile=profile, figures=figures)

    @staticmethod
    def _load_workflow(
        experiment: Experiment,
    ) -> Callable[..., None]:
        module_name = f"_satute_{experiment.id}_workflow"
        specification = importlib.util.spec_from_file_location(
            module_name, experiment.workflow_path
        )
        if specification is None or specification.loader is None:
            raise ExperimentError(f"Cannot load workflow for {experiment.id}")
        module = importlib.util.module_from_spec(specification)
        specification.loader.exec_module(module)
        if not callable(getattr(module, "run", None)):
            raise ExperimentError(
                f"{experiment.workflow_path} must define run(context, *, profile, figures)"
            )
        return module.run


__all__ = [
    "Experiment",
    "ExperimentCatalog",
    "ExperimentError",
    "WorkflowContext",
]
