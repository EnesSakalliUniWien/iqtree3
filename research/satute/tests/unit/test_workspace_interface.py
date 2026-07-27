#!/usr/bin/env python3
"""Behavioral contract for the single workspace interface and environment."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[2]
ENTRY_POINT = PROJECT / "satute"
EXPECTED_PYTHON = PROJECT / ".venv" / "bin" / "python"


def main() -> None:
    if not ENTRY_POINT.is_file():
        raise SystemExit(f"Missing workspace entry point: {ENTRY_POINT}")
    if not EXPECTED_PYTHON.is_file():
        raise SystemExit(f"Missing canonical virtual environment: {EXPECTED_PYTHON}")

    result = subprocess.run(
        [str(ENTRY_POINT), "doctor", "--json"],
        cwd=PROJECT,
        check=True,
        capture_output=True,
        text=True,
    )
    report = json.loads(result.stdout)
    if Path(report["entry_point"]).resolve() != ENTRY_POINT.resolve():
        raise SystemExit(f"Doctor reported the wrong entry point: {report}")
    if Path(report["python"]).resolve() != EXPECTED_PYTHON.resolve():
        raise SystemExit(f"Doctor escaped the canonical virtual environment: {report}")
    if report["virtual_environments"] != [str(PROJECT / ".venv")]:
        raise SystemExit(f"Workspace has multiple virtual environments: {report}")

    help_result = subprocess.run(
        [str(ENTRY_POINT), "--help"],
        cwd=PROJECT,
        check=True,
        capture_output=True,
        text=True,
    )
    for command in (
        "bootstrap",
        "doctor",
        "data",
        "experiment",
        "cluster",
        "publication-figures",
        "check",
        "manuscript",
        "developer",
    ):
        if command not in help_result.stdout:
            raise SystemExit(f"Workspace interface omits command: {command}")

    makefile = (PROJECT / "Makefile").read_text(encoding="utf-8")
    for forbidden in ("miniconda", "/usr/bin/python", "command -v python3"):
        if forbidden in makefile:
            raise SystemExit(f"Makefile still discovers another interpreter: {forbidden}")
    alternate_dependency_files = (
        PROJECT / "requirements.txt",
        PROJECT / "Pipfile",
        PROJECT / "environment.yml",
        PROJECT / "poetry.lock",
    )
    existing_alternates = [path.name for path in alternate_dependency_files if path.exists()]
    if existing_alternates:
        raise SystemExit(
            f"Dependency truth is duplicated outside pyproject.toml: {existing_alternates}"
        )

    active_guides = tuple(
        path
        for path in PROJECT.rglob("*.md")
        if not any(
            part in {"artifacts", "archive", ".venv", "tmp"}
            for part in path.relative_to(PROJECT).parts
        )
        and "tools/cluster" not in path.relative_to(PROJECT).as_posix()
    )
    for guide in active_guides:
        text = guide.read_text(encoding="utf-8")
        for forbidden in (
            "python generate_experiment_data.py",
            "python ../../generate_experiment_data.py",
            "make status",
            "make smoke",
            "make check",
            "make experiment",
            "make publication-figures",
            "make manuscript",
            "make audit",
            "make -C",
            "python3 research/satute/",
            "${PYTHON_BIN:-python3}",
        ):
            if forbidden in text:
                raise SystemExit(
                    f"Active guide bypasses ./satute: {guide.relative_to(PROJECT)}: "
                    f"{forbidden}"
                )

    print("Single workspace interface and virtual environment: passed")


if __name__ == "__main__":
    main()
