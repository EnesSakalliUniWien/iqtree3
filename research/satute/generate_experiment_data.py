#!/usr/bin/env python3
"""Generate every managed SatuTe experiment dataset through one entry point."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


PROJECT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT / "src"))

from satute_analysis.experiments import ExperimentCatalog, ExperimentError


def selected_runs(catalog, requested):
    available = {
        f"{experiment.id}/{study.id}": (experiment, study)
        for experiment in catalog
        for study in experiment.studies
        if study.scale != "manual"
    }
    if not requested:
        return list(available.values())

    unknown = sorted(set(requested) - available.keys())
    if unknown:
        choices = "\n  ".join(available)
        raise ExperimentError(
            f"Unknown experiment/study selection: {unknown}\nAvailable:\n  {choices}"
        )
    return [available[name] for name in requested]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--profile",
        choices=("smoke", "full"),
        help="smoke validates wiring; full generates the declared evidence grids.",
    )
    parser.add_argument(
        "--map",
        action="store_true",
        help="Print the complete simulation, analysis, and artifact-lineage map.",
    )
    parser.add_argument(
        "--only",
        action="append",
        metavar="EXPERIMENT/STUDY",
        help="Run only this study; repeat for multiple selective reruns.",
    )
    parser.add_argument(
        "--plan",
        action="store_true",
        help="Print the complete run plan without executing drivers.",
    )
    parser.add_argument(
        "--skip-rate-category-check",
        action="store_true",
        help="Skip the native +G4/+R4/+I rate-category preflight.",
    )
    parser.add_argument(
        "--figures",
        action="store_true",
        help="After data generation, build all publication power figures from the canonical case map.",
    )
    parser.add_argument(
        "--figure-mode",
        choices=("render", "index-existing"),
        default="render",
        help="Render figures from source tables, or organize already approved manuscript PDFs.",
    )
    args = parser.parse_args()

    if args.map:
        map_path = PROJECT / "docs/overview/simulation-analysis-map.md"
        print(map_path.read_text(encoding="utf-8"), end="")
        return 0
    if args.profile is None and not args.figures:
        parser.error("--profile is required unless --map or --figures is used")

    catalog = None
    runs = []
    if args.profile is not None:
        try:
            catalog = ExperimentCatalog.discover(PROJECT)
            runs = selected_runs(catalog, args.only)
        except ExperimentError as error:
            parser.error(str(error))

    if args.profile is not None:
        print(f"SatuTe data plan: {args.profile}")
    if args.profile is not None and not args.skip_rate_category_check:
        print(
            "  00. native/rate_categories [correctness gate; local]\n"
            "      Verify IQ-TREE category assignments, multipliers, effective "
            "lengths, and pooled statistics."
        )
    for index, (experiment, study) in enumerate(runs, start=1):
        print(
            f"  {index:02d}. {experiment.id}/{study.id} "
            f"[{study.outcome}; {study.scale}]"
        )
        print(f"      {study.question}")
    manual = [
        f"{experiment.id}/{study.id}"
        for experiment in (catalog or ())
        for study in experiment.studies
        if study.scale == "manual"
    ]
    if manual:
        print("  Manual inputs required: " + ", ".join(manual))

    if args.plan:
        if args.figures:
            subprocess.run(
                [
                    sys.executable,
                    str(PROJECT / "tools/maintenance/generate_publication_power_figures.py"),
                    "--mode",
                    args.figure_mode,
                    "--plan",
                ],
                cwd=PROJECT,
                check=True,
            )
        return 0

    if args.profile is not None and not args.skip_rate_category_check:
        iqtree = PROJECT.parents[1] / "build/iqtree3"
        rate_output = PROJECT / "artifacts/local/check/rate-categories"
        print("\n[preflight] Verifying native rate-category handling")
        subprocess.run(
            [
                sys.executable,
                str(PROJECT / "tests/native/test_satute_rate_categories.py"),
                str(iqtree),
                str(rate_output),
            ],
            cwd=PROJECT,
            check=True,
        )

    for index, (experiment, study) in enumerate(runs, start=1):
        print(f"\n[{index}/{len(runs)}] Generating {experiment.id}/{study.id}")
        try:
            catalog.run(
                experiment.id,
                study=study.id,
                profile=args.profile,
                python=sys.executable,
            )
        except ExperimentError as error:
            parser.error(str(error))
    if args.figures:
        subprocess.run(
            [
                sys.executable,
                str(PROJECT / "tools/maintenance/generate_publication_power_figures.py"),
                "--mode",
                args.figure_mode,
            ],
            cwd=PROJECT,
            check=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
