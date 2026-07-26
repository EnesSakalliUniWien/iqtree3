#!/usr/bin/env python3
"""Single entry point: reproduce every scientific result this machine can
actually produce, and say plainly what it cannot.

Of the 9 registered experiment studies, only 3 are scale="local":
  001/historical_power_reconstruction, 003/invariant_mixture_power,
  004/saturation_diagnostics
5 are scale="cluster" (all of 002 -- gtr_power, protein_power,
misspecification_power, fdr_exact_null, fdr_mixed_null), and 005 is
scale="manual" (needs explicit biological alignment inputs). The
misspecification figure contract (MISSPECIFICATION_CONTRACT in
publication_figures.py) requires both tree cases x all 3 evaluation models x
4 series -- a real cluster-scale grid, not something to fake at reduced
replicate count locally. Running this script does NOT run those five
studies; it says so, and tells you the command that does
(``./satute cluster submit``).

What it DOES run, in order:
  1. doctor                          -- environment sanity
  2. native rate-category preflight  -- correctness gate IQ-TREE-side
  3. the 3 local-scale registered studies, at --profile (default full)
  4. validate_weighted_formula.py's full experiment suite -- the actual
     scientific payoff of the misspecification investigation: calibration,
     power, and misspecification-robustness of the (now centred)
     eigenvalue_weighted statistic. This is genuinely local-scale and is
     the direct local analogue of 002's misspecification_power study.

Usage:
    research/satute/satute developer reproduce-scientific-results -- \\
        --iqtree ./build/iqtree3

    # fast iteration, skip the slow bits:
    ... --profile smoke --reps-scale 0.1
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
REPOSITORY = PROJECT.parents[1]
sys.path.insert(0, str(PROJECT / "src"))

from satute_analysis.experiments import ExperimentCatalog, ExperimentError  # noqa: E402


def section(title: str) -> None:
    print(f"\n{'#' * 70}\n# {title}\n{'#' * 70}")


def stage_doctor() -> dict:
    section("1. Environment")
    result = subprocess.run(["./satute", "doctor"], cwd=PROJECT, capture_output=True, text=True)
    print(result.stdout, end="")
    ok = result.returncode == 0
    print("PASS: workspace environment healthy" if ok else "FAIL: workspace environment inconsistent -- run ./satute bootstrap")
    return {"ok": ok}


def stage_rate_category_preflight(iqtree: str) -> dict:
    section("2. Native rate-category preflight")
    outdir = PROJECT / "artifacts/local/check/rate-categories"
    started = time.time()
    result = subprocess.run(
        [sys.executable, str(PROJECT / "tests/native/test_satute_rate_categories.py"), iqtree, str(outdir)],
        cwd=PROJECT, check=False,
    )
    ok = result.returncode == 0
    elapsed = time.time() - started
    print(f"{'PASS' if ok else 'FAIL'} ({elapsed:.1f}s)")
    return {"ok": ok, "elapsed_s": elapsed}


def stage_local_experiments(profile: str) -> dict:
    section("3. Local-scale registered experiments")
    catalog = ExperimentCatalog.discover(PROJECT)
    local_studies = [
        (experiment, study)
        for experiment in catalog
        for study in experiment.studies
        if study.scale == "local"
    ]
    cluster_studies = [
        (experiment, study)
        for experiment in catalog
        for study in experiment.studies
        if study.scale == "cluster"
    ]
    manual_studies = [
        (experiment, study)
        for experiment in catalog
        for study in experiment.studies
        if study.scale == "manual"
    ]

    print(f"Running {len(local_studies)} local-scale studies at --profile {profile}:")
    results = {}
    for experiment, study in local_studies:
        key = f"{experiment.id}/{study.id}"
        print(f"\n  [{key}] {study.question}")
        started = time.time()
        try:
            catalog.run(experiment.id, study=study.id, profile=profile, python=sys.executable)
            results[key] = {"ok": True, "elapsed_s": time.time() - started}
            print(f"  PASS ({results[key]['elapsed_s']:.1f}s)")
        except ExperimentError as error:
            results[key] = {"ok": False, "error": str(error), "elapsed_s": time.time() - started}
            print(f"  FAIL: {error}")

    if cluster_studies:
        print(f"\nNOT run here -- {len(cluster_studies)} cluster-scale studies "
              "(need the LiSC campaign, not this laptop):")
        for experiment, study in cluster_studies:
            print(f"  {experiment.id}/{study.id} -- {study.question}")
        print("  Run: ./satute cluster plan   (then)   ./satute cluster submit")

    if manual_studies:
        print(f"\nNOT run here -- {len(manual_studies)} manual-input studies:")
        for experiment, study in manual_studies:
            print(f"  {experiment.id}/{study.id} -- needs explicit inputs, see "
                  f"{experiment.directory / 'README.md'}")

    return {"studies": results, "skipped_cluster": len(cluster_studies), "skipped_manual": len(manual_studies)}


def stage_weighted_formula_validation(iqtree: str, reps_scale: float, which: str) -> dict:
    section("4. Weighted-formula scientific validation "
            "(calibration, power, misspecification-robustness)")
    print("Local analogue of 002/misspecification_power: same true/misspecified model\n"
          "design (GTR-true vs K2P-misspecified), run at local scale with self-checking\n"
          "assertions rather than the cluster-scale publication grid.\n")
    outdir = PROJECT / "artifacts/local/scientific-results/validate-weighted-formula"
    started = time.time()
    result = subprocess.run(
        [sys.executable, str(PROJECT / "tests/native/validate_weighted_formula.py"),
         "--iqtree", iqtree, "--outdir", str(outdir), "--which", which, "--reps-scale", str(reps_scale)],
        cwd=PROJECT, check=False,
    )
    ok = result.returncode == 0
    elapsed = time.time() - started
    print(f"\n{'PASS' if ok else 'FAIL'} ({elapsed:.1f}s) -- outputs in {outdir}")
    return {"ok": ok, "elapsed_s": elapsed, "outdir": str(outdir)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--iqtree", default=str(REPOSITORY / "build/iqtree3"))
    parser.add_argument("--profile", choices=("smoke", "full"), default="full",
                         help="Profile for the 3 local-scale registered experiments.")
    parser.add_argument("--reps-scale", type=float, default=1.0,
                         help="Passed through to validate_weighted_formula.py.")
    parser.add_argument("--which", default="all",
                         help="Passed through to validate_weighted_formula.py's --which.")
    parser.add_argument("--skip", action="append", default=[],
                         choices=("doctor", "preflight", "experiments", "validation"),
                         help="Skip a stage; repeatable.")
    args = parser.parse_args()

    iqtree = str(Path(args.iqtree))
    started = time.time()
    report: dict = {"stages": {}}

    if "doctor" not in args.skip:
        report["stages"]["doctor"] = stage_doctor()
    if "preflight" not in args.skip:
        report["stages"]["preflight"] = stage_rate_category_preflight(iqtree)
    if "experiments" not in args.skip:
        report["stages"]["local_experiments"] = stage_local_experiments(args.profile)
    if "validation" not in args.skip:
        report["stages"]["weighted_formula_validation"] = stage_weighted_formula_validation(
            iqtree, args.reps_scale, args.which
        )

    section("Summary")
    failures = []
    for name, stage in report["stages"].items():
        if name == "local_experiments":
            for key, study_result in stage.get("studies", {}).items():
                if not study_result.get("ok"):
                    failures.append(f"local_experiments/{key}")
        elif not stage.get("ok", True):
            failures.append(name)

    elapsed_total = time.time() - started
    report["elapsed_s"] = elapsed_total
    report["failures"] = failures

    summary_path = PROJECT / "artifacts/local/scientific-results/summary.json"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(report, indent=2, default=str))

    print(f"Total time: {elapsed_total:.1f}s")
    print(f"Summary written to {summary_path}")
    if failures:
        print(f"\n{len(failures)} stage(s) FAILED: {failures}")
        return 1
    print("\nAll runnable stages passed. Cluster-scale studies still need "
          "./satute cluster submit for publication-grade numbers.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
