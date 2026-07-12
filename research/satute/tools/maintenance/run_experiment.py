#!/usr/bin/env python3
"""Registry-aware experiment entry point."""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[2]
REPOSITORY = PROJECT.parents[1]
IQTREE = REPOSITORY / "build/iqtree3"


def run(command):
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(PROJECT / "src") + os.pathsep + environment.get("PYTHONPATH", "")
    print("+", " ".join(map(str, command)))
    subprocess.run([str(item) for item in command], cwd=PROJECT, env=environment, check=True)


def paper_reconstruction(profile, figures):
    output = PROJECT / "artifacts/local/paper-reconstruction"
    if profile == "smoke":
        output = PROJECT / "artifacts/local/smoke/001_paper_reconstruction"
    if figures:
        run(
            [
                sys.executable,
                "experiments/001_paper_reconstruction/render.py",
                "--summary",
                output / "fig2_reconstruction_summary.tsv",
                "--output",
                output / "figures/fig2_reconstruction.svg",
            ]
        )
        return
    command = [
        sys.executable,
        "experiments/001_paper_reconstruction/run.py",
        "--iqtree",
        IQTREE,
        "--outdir",
        output,
    ]
    if profile == "smoke":
        command.extend(
            ["--reps", "1", "--site-lengths", "100", "--branch-lengths", "0.5,1.0", "--tree-cases", "five_external,sixteen_internal"]
        )
    else:
        command.extend(["--reps", "1000", "--paper-grid"])
    run(command)


def relative_weighting(profile, figures):
    output = PROJECT / "artifacts/local/relative-weighting"
    if profile == "smoke":
        output = PROJECT / "artifacts/local/smoke/002_relative_weighting"
    if figures:
        command = [
            sys.executable,
            "experiments/002_relative_weighting/render.py",
            "--summary",
            output / "head_to_head_summary.tsv",
            "--outdir",
            output / "figures",
            "--expected-reps",
            "1" if profile == "smoke" else "1000",
            "--tree-cases",
            "five_external,sixteen_internal",
        ]
        if profile == "smoke":
            command.append("--allow-incomplete")
        run(command)
        return
    command = [
        sys.executable,
        "experiments/002_relative_weighting/run.py",
        "--iqtree",
        IQTREE,
        "--outdir",
        output,
        "--simulation-models",
        "JC",
        "--evaluation-models",
        "JC",
    ]
    if profile == "smoke":
        command.extend(
            [
                "--reps",
                "1",
                "--site-lengths",
                "100",
                "--branch-lengths",
                "0.5,1.0",
                "--tree-cases",
                "five_external,sixteen_internal",
                "--evonaps-branch-lengths",
                PROJECT / "references/evonaps/evonaps_16taxon_branch_lengths.tsv",
            ]
        )
    else:
        command.extend(
            [
                "--reps",
                "1000",
                "--paper-grid",
                "--tree-cases",
                "five_external,sixteen_internal",
                "--evonaps-branch-lengths",
                PROJECT / "references/evonaps/evonaps_16taxon_branch_lengths.tsv",
            ]
        )
    run(command)


def invariant(profile, figures):
    output = PROJECT / "artifacts/local/weighted-invariant-full"
    if profile == "smoke":
        output = PROJECT / "artifacts/local/smoke/003_invariant_mixture"
    if figures:
        run(
            [
                sys.executable,
                "experiments/003_invariant_mixture/render_full.py",
                "--input-dir",
                output,
                "--output-prefix",
                output / "figures/full-analysis",
                "--png-scale",
                "3",
            ]
        )
        return
    command = [
        sys.executable,
        "experiments/003_invariant_mixture/run.py",
        "--outdir",
        output,
        "--figure-prefix",
        output / "figures/overview",
        "--png-scale",
        "3",
    ]
    if profile == "smoke":
        command.extend(
            [
                "--site-lengths",
                "100,1000",
                "--branch-lengths",
                "20,50,10000",
                "--null-reps",
                "200",
                "--alternative-reps",
                "200",
            ]
        )
    else:
        command.extend(["--null-reps", "20000", "--alternative-reps", "20000"])
    run(command)


def saturation_diagnostics(profile, figures):
    output = PROJECT / "artifacts/local/saturation-diagnostics"
    if profile == "smoke":
        output = PROJECT / "artifacts/local/smoke/004_saturation_diagnostics"
    command = [
        sys.executable,
        "experiments/004_saturation_diagnostics/run.py",
        "--iqtree",
        IQTREE,
        "--outdir",
        output,
        "--figure-prefix",
        output / "figures/main",
        "--category-figure-prefix",
        output / "figures/categories",
    ]
    if profile == "smoke":
        command.extend(
            [
                "--sites",
                "100",
                "--branch-lengths",
                "0.5,10",
                "--null-length",
                "10",
                "--reps",
                "2",
                "--null-reps",
                "10",
                "--workers",
                "2",
            ]
        )
    if figures:
        command.append("--plot-only")
    run(command)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", required=True)
    parser.add_argument("--profile", choices=("smoke", "full"), default="smoke")
    parser.add_argument("--figures", action="store_true")
    args = parser.parse_args()
    registry = json.loads((PROJECT / "experiments/registry.yaml").read_text())
    known = {entry["id"] for entry in registry["experiments"]}
    if args.experiment not in known:
        raise SystemExit(f"Unknown experiment {args.experiment}; choose from {sorted(known)}")
    dispatch = {
        "001_paper_reconstruction": paper_reconstruction,
        "002_relative_weighting": relative_weighting,
        "003_invariant_mixture": invariant,
        "004_saturation_diagnostics": saturation_diagnostics,
    }
    if args.experiment in dispatch:
        dispatch[args.experiment](args.profile, args.figures)
        return
    raise SystemExit(
        "005_biological_windows requires alignment inputs and is intentionally not run "
        "without them; see experiments/005_biological_windows/README.md"
    )


if __name__ == "__main__":
    main()
