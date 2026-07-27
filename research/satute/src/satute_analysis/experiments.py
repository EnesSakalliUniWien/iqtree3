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
STUDY_ID = re.compile(r"^[a-z][a-z0-9]*(?:_[a-z0-9]+)*$")
REQUIRED_FIELDS = {
    "schema_version",
    "id",
    "question",
    "status",
    "driver",
    "consumers",
    "studies",
    "default_study",
    "config",
}
OUTPUT_SCOPES = ("artifacts/local", "artifacts/cluster")
RELEASE_SCOPE = "artifacts/releases"


class ExperimentError(RuntimeError):
    """Raised when an experiment definition or workflow is invalid."""


def _contained_path(root: Path, relative_path: str | Path, label: str) -> Path:
    """Resolve a relative path and reject escapes from its owning directory."""

    path = Path(relative_path)
    if path.is_absolute():
        raise ExperimentError(f"{label} must be relative to {root}: {relative_path}")
    resolved = (root / path).resolve()
    try:
        resolved.relative_to(root.resolve())
    except ValueError as error:
        raise ExperimentError(f"{label} escapes {root}: {relative_path}") from error
    return resolved


def _dedicated_artifact_path(
    project_root: Path,
    relative_path: str,
    *,
    label: str,
    scopes: tuple[str, ...],
) -> Path:
    """Require a dedicated child directory in one of the declared artifact scopes."""

    resolved = _contained_path(project_root, relative_path, label)
    allowed_roots = tuple((project_root / scope).resolve() for scope in scopes)
    for allowed_root in allowed_roots:
        try:
            child = resolved.relative_to(allowed_root)
        except ValueError:
            continue
        if child.parts:
            return resolved
    scope_list = ", ".join(scopes)
    raise ExperimentError(
        f"{label} must be a dedicated directory under one of: {scope_list}"
    )


@dataclass(frozen=True)
class ExperimentStudy:
    """One named scientific study exposed by an experiment module."""

    id: str
    question: str
    outcome: str
    scale: str


@dataclass(frozen=True)
class Experiment:
    """Validated, immutable description of one experiment module."""

    directory: Path
    id: str
    question: str
    status: str
    consumers: tuple[str, ...]
    studies: tuple[ExperimentStudy, ...]
    default_study_id: str
    driver_path: Path
    workflow_path: Path | None
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

    @property
    def study_ids(self) -> tuple[str, ...]:
        return tuple(study.id for study in self.studies)

    def study(self, study_id: str | None = None) -> ExperimentStudy:
        selected_id = study_id or self.default_study_id
        for study in self.studies:
            if study.id == selected_id:
                return study
        raise ExperimentError(
            f"Unknown study {selected_id!r} for {self.id}; "
            f"choose from {list(self.study_ids)}"
        )

    def file(self, relative_path: str | Path) -> Path:
        """Resolve a file owned by this experiment without escaping its module."""

        if Path(relative_path).is_absolute():
            raise ExperimentError(
                f"{self.id} path must be relative to its experiment module: {relative_path}"
            )
        return _contained_path(
            self.directory,
            relative_path,
            f"{self.id} path",
        )


@dataclass(frozen=True)
class _WorkflowContext:
    """Capabilities supplied to experiment workflow adapters."""

    project_root: Path
    experiment: Experiment
    study: ExperimentStudy
    python: Path

    @property
    def repository_root(self) -> Path:
        return self.project_root.parents[1]

    @property
    def iqtree(self) -> Path:
        return self.repository_root / "build/iqtree3"

    def output_dir(self, profile: str) -> Path:
        study_suffix = (
            Path(self.study.id) if len(self.experiment.studies) > 1 else Path()
        )
        if profile == "smoke":
            return (
                self.project_root
                / "artifacts/local/smoke"
                / self.experiment.id
                / study_suffix
                / f"schema-{self.experiment.config_schema_version}"
            )
        return self.experiment.output_path / study_suffix

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
        numbered_directories = (
            {
                directory.resolve()
                for directory in experiment_root.iterdir()
                if directory.is_dir() and EXPERIMENT_ID.fullmatch(directory.name)
            }
            if experiment_root.is_dir()
            else set()
        )
        manifest_paths = sorted(experiment_root.glob(f"*/{MANIFEST_NAME}"))
        manifest_directories = {path.parent.resolve() for path in manifest_paths}
        missing_manifests = sorted(numbered_directories - manifest_directories)
        if missing_manifests:
            raise ExperimentError(
                "Numbered experiment directories are missing experiment.json: "
                + ", ".join(path.name for path in missing_manifests)
            )
        experiments = tuple(cls._read_manifest(path) for path in manifest_paths)
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
        execution = manifest.get("execution", "managed")
        if execution not in {"managed", "manual"}:
            raise ExperimentError(
                f"{experiment_id} execution must be either 'managed' or 'manual'"
            )
        workflow_value = manifest.get("workflow")
        if execution == "managed":
            workflow = cls._required_text(manifest, "workflow", path)
        elif workflow_value is not None:
            raise ExperimentError(
                f"{experiment_id} manual experiments must omit the workflow field"
            )
        else:
            workflow = None

        if not isinstance(manifest["consumers"], list):
            raise ExperimentError(f"{experiment_id} consumers must be a list")
        consumers = tuple(manifest["consumers"])
        if any(
            not isinstance(consumer, str) or not consumer.strip()
            for consumer in consumers
        ):
            raise ExperimentError(f"{experiment_id} consumers must contain non-empty strings")
        if len(consumers) != len(set(consumers)):
            raise ExperimentError(f"{experiment_id} consumers must be unique")

        studies_value = manifest["studies"]
        if not isinstance(studies_value, list) or not studies_value:
            raise ExperimentError(f"{experiment_id} studies must be a non-empty list")
        studies = tuple(
            cls._read_study(experiment_id, value, path) for value in studies_value
        )
        study_ids = tuple(study.id for study in studies)
        if len(study_ids) != len(set(study_ids)):
            raise ExperimentError(f"{experiment_id} study ids must be unique")
        default_study_id = cls._required_text(manifest, "default_study", path)
        if default_study_id not in study_ids:
            raise ExperimentError(
                f"{experiment_id} default_study must be one of {list(study_ids)}"
            )

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
        output_path = _dedicated_artifact_path(
            project_root,
            output,
            label=f"{experiment_id} output directory",
            scopes=OUTPUT_SCOPES,
        )

        release_value = manifest.get("release")
        if release_value is not None and (
            not isinstance(release_value, str) or not release_value.strip()
        ):
            raise ExperimentError(f"{experiment_id} release must be a non-empty string")
        release_path = (
            _dedicated_artifact_path(
                project_root,
                release_value,
                label=f"{experiment_id} release",
                scopes=(RELEASE_SCOPE,),
            )
            if release_value
            else None
        )
        directory = path.parent.resolve()
        readme_path = directory / "README.md"
        if not readme_path.is_file():
            raise ExperimentError(f"{experiment_id} README is missing: {readme_path}")

        experiment = Experiment(
            directory=directory,
            id=experiment_id,
            question=question,
            status=status,
            consumers=consumers,
            studies=studies,
            default_study_id=default_study_id,
            driver_path=_contained_path(
                directory, driver, f"{experiment_id} driver path"
            ),
            workflow_path=(
                _contained_path(
                    directory, workflow, f"{experiment_id} workflow path"
                )
                if workflow
                else None
            ),
            output_path=output_path,
            config_schema_version=config_schema_version,
            release_path=release_path,
            _configuration=deepcopy(configuration),
        )
        owned_files = [("driver", experiment.driver_path)]
        if experiment.workflow_path is not None:
            owned_files.append(("workflow", experiment.workflow_path))
        for role, owned_path in owned_files:
            if not owned_path.is_file():
                raise ExperimentError(f"{experiment.id} {role} is missing: {owned_path}")
        return experiment

    @classmethod
    def _read_study(
        cls,
        experiment_id: str,
        value: Any,
        source: Path,
    ) -> ExperimentStudy:
        if not isinstance(value, dict):
            raise ExperimentError(f"{experiment_id} study entries must be objects")
        study_id = cls._required_text(value, "id", source)
        if not STUDY_ID.fullmatch(study_id):
            raise ExperimentError(f"{experiment_id} has invalid study id: {study_id}")
        question = cls._required_text(value, "question", source)
        outcome = cls._required_text(value, "outcome", source)
        if outcome not in {
            "baseline",
            "power",
            "calibration",
            "diagnostic",
            "application",
        }:
            raise ExperimentError(
                f"{experiment_id}/{study_id} has invalid outcome: {outcome}"
            )
        scale = cls._required_text(value, "scale", source)
        if scale not in {"local", "cluster", "manual"}:
            raise ExperimentError(
                f"{experiment_id}/{study_id} has invalid scale: {scale}"
            )
        return ExperimentStudy(study_id, question, outcome, scale)

    @staticmethod
    def _required_text(
        values: Mapping[str, Any], field_name: str, source: Path
    ) -> str:
        value = values.get(field_name)
        if not isinstance(value, str) or not value.strip():
            raise ExperimentError(
                f"{source} field {field_name!r} must be a non-empty string"
            )
        return value

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
        study: str | None = None,
        python: Path | str = sys.executable,
    ) -> None:
        if profile not in {"smoke", "full"}:
            raise ExperimentError(f"Unknown profile {profile!r}; choose smoke or full")
        experiment = self.get(experiment_id)
        selected_study = experiment.study(study)
        if experiment.workflow_path is None:
            raise ExperimentError(
                f"{experiment.id} requires explicit inputs and cannot use the common "
                f"smoke/full workflow; see {experiment.directory / 'README.md'}"
            )
        workflow = self._load_workflow(experiment)
        context = _WorkflowContext(
            self.project_root,
            experiment,
            selected_study,
            Path(python),
        )
        workflow(
            context,
            profile=profile,
            figures=figures,
            study=selected_study.id,
        )

    @staticmethod
    def _load_workflow(
        experiment: Experiment,
    ) -> Callable[..., None]:
        if experiment.workflow_path is None:
            raise ExperimentError(f"{experiment.id} has no managed workflow")
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
                f"{experiment.workflow_path} must define "
                "run(context, *, profile, figures, study)"
            )
        return module.run


__all__ = [
    "Experiment",
    "ExperimentCatalog",
    "ExperimentError",
    "ExperimentStudy",
]
