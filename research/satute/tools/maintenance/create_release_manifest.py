#!/usr/bin/env python3
"""Create or refresh checksummed release provenance."""

import argparse
import json
import sys
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT / "src"))

from satute_analysis.artifacts import write_manifest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", required=True)
    parser.add_argument("--version", required=True)
    args = parser.parse_args()
    release_name = args.experiment.removeprefix("003_").replace("_", "-") + "-" + args.version
    release = PROJECT / "artifacts/releases" / release_name
    if not release.exists():
        raise SystemExit(f"Release directory does not exist: {release}")
    config_path = PROJECT / "config/experiments" / f"{args.experiment}.yaml"
    config = json.loads(config_path.read_text())
    outputs = sorted(path for path in release.rglob("*") if path.is_file() and path.name != "manifest.json")
    write_manifest(
        release / "manifest.json",
        args.experiment,
        config,
        PROJECT.parents[1],
        outputs,
    )
    print(release / "manifest.json")


if __name__ == "__main__":
    main()
