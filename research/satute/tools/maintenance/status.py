#!/usr/bin/env python3
"""Print one compact, observable view of the SatuTe workspace."""

import json
import subprocess
from datetime import datetime
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[2]
REPOSITORY = PROJECT.parents[1]


def directory_size(path):
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def human_size(value):
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if value < 1024.0 or unit == "TB":
            return f"{value:.1f} {unit}"
        value /= 1024.0


def git_lines():
    result = subprocess.run(
        ["git", "status", "--short"],
        cwd=REPOSITORY,
        text=True,
        stdout=subprocess.PIPE,
        check=True,
    )
    return [line for line in result.stdout.splitlines() if line]


def main():
    registry = json.loads((PROJECT / "experiments/registry.yaml").read_text())
    print("SatuTe workspace")
    print(f"  root: {PROJECT}")
    print("\nExperiments")
    for experiment in registry["experiments"]:
        release = experiment.get("release", "-")
        print(f"  {experiment['id']:<28} {experiment['status']:<20} {release}")

    print("\nArtifacts")
    for name in ("local", "cluster", "releases"):
        path = PROJECT / "artifacts" / name
        print(f"  {name:<10} {human_size(directory_size(path))}")

    manuscript = PROJECT / "manuscript" / "manuscript.pdf"
    print("\nManuscript")
    if manuscript.exists():
        stamp = datetime.fromtimestamp(manuscript.stat().st_mtime).isoformat(timespec="seconds")
        print(f"  built: {stamp}; {human_size(manuscript.stat().st_size)}")
    else:
        print("  not built")

    changes = git_lines()
    production = [line for line in changes if "tree/satute.cpp" in line]
    research = [line for line in changes if "research/satute" in line or "doc/satute" in line]
    other = [line for line in changes if line not in production and line not in research]
    print("\nGit changes")
    print(f"  production: {len(production)}")
    print(f"  research:   {len(research)}")
    print(f"  other:      {len(other)}")
    for line in production + other[:8]:
        print(f"    {line}")


if __name__ == "__main__":
    main()
