#!/usr/bin/env python3
"""Audit experiment modules, project seams, and manuscript dependencies."""

import json
import os
import re
import subprocess
import sys
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[2]
REPOSITORY = PROJECT.parents[1]
sys.path.insert(0, str(PROJECT / "src"))

from satute_analysis.artifacts import file_sha256
from satute_analysis.experiments import ExperimentCatalog, ExperimentError


def fail(message, failures):
    failures.append(message)
    print(f"FAIL: {message}")


def main():
    failures = []
    try:
        catalog = ExperimentCatalog.discover(PROJECT)
    except ExperimentError as error:
        fail(f"experiment discovery failed: {error}", failures)
        catalog = ()
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(PROJECT / "src") + os.pathsep + environment.get(
        "PYTHONPATH", ""
    )
    for experiment in catalog:
        if experiment.driver_path.exists():
            result = subprocess.run(
                [sys.executable, str(experiment.driver_path), "--help"],
                cwd=PROJECT,
                env=environment,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                text=True,
            )
            if result.returncode:
                fail(
                    f"driver CLI failed: {experiment.id}: {result.stderr[-500:]}",
                    failures,
                )

    workspace_interface = PROJECT / "satute"
    result = subprocess.run(
        [str(workspace_interface), "data", "--profile", "smoke", "--plan"],
        cwd=PROJECT,
        env=environment,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode:
        fail(f"workspace data-generation plan failed: {result.stderr[-500:]}", failures)

    simulation_map_path = PROJECT / "docs/overview/simulation-analysis-map.md"
    simulation_map = simulation_map_path.read_text(encoding="utf-8")
    required_map_entries = (
        "experiments/001_paper_reconstruction/run.py",
        "experiments/002_relative_weighting/run.py",
        "experiments/002_relative_weighting/fdr_calibration_run.py",
        "experiments/003_invariant_mixture/run.py",
        "experiments/004_saturation_diagnostics/run.py",
        "experiments/005_biological_windows/analyze_dataset.py",
        "experiments/005_biological_windows/prepare_hcv_example.py",
        "tests/exact/test_saturation_null.py",
        "tests/native/compare_formulas.py",
        "tests/native/test_satute_rate_categories.py",
        "tests/native/test_satute_reference_general.py",
        "tests/native/test_satute_simulated.sh",
        "src/satute_analysis/statistics.py",
        "experiments/002_relative_weighting/combine.py",
        "experiments/002_relative_weighting/fdr_calibration_postprocess.py",
        "experiments/002_relative_weighting/recalibrate_heldout_from_truth.py",
        "experiments/002_relative_weighting/render_power_summary.py",
        "experiments/003_invariant_mixture/render_full.py",
        "experiments/004_saturation_diagnostics/summarize_null.py",
        "experiments/005_biological_windows/render_hcv_example.py",
    )
    missing_map_entries = [
        entry for entry in required_map_entries if entry not in simulation_map
    ]
    if missing_map_entries:
        fail(
            "simulation/analysis map omits internal adapters: "
            + ", ".join(missing_map_entries),
            failures,
        )

    manuscript = PROJECT / "manuscript"
    tex_files = [manuscript / "main.tex", *sorted((manuscript / "sections").rglob("*.tex"))]
    for tex_file in tex_files:
        text = tex_file.read_text(encoding="utf-8")
        for match in re.finditer(r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}", text):
            asset = manuscript / match.group(1)
            if not asset.exists():
                fail(f"missing manuscript asset: {asset}", failures)
            ignored = subprocess.run(
                ["git", "check-ignore", "-q", str(asset)], cwd=REPOSITORY
            ).returncode == 0
            if ignored:
                fail(f"manuscript asset is ignored: {asset}", failures)

    peer_imports = []
    for path in [*PROJECT.glob("experiments/**/*.py"), *PROJECT.glob("tests/**/*.py")]:
        text = path.read_text(encoding="utf-8")
        if "head_to_head_satute_simulations.py" in text:
            peer_imports.append(path)
    if peer_imports:
        fail("peer-script imports remain: " + ", ".join(map(str, peer_imports)), failures)

    for script in sorted(PROJECT.glob("tools/cluster/**/*")):
        if script.suffix not in {".sh", ".slurm"}:
            continue
        result = subprocess.run(
            ["bash", "-n", str(script)], stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True
        )
        if result.returncode:
            fail(f"shell syntax failed: {script}: {result.stderr.strip()}", failures)

    compatibility = {
        REPOSITORY / "doc/satute-wiki": PROJECT,
        REPOSITORY / "doc/satute-phase1": PROJECT / "tests/native",
    }
    for link, expected in compatibility.items():
        if not link.is_symlink() or link.resolve() != expected.resolve():
            fail(f"compatibility link invalid: {link} -> {expected}", failures)

    for manifest in sorted((PROJECT / "artifacts/releases").glob("*/manifest.json")):
        payload = json.loads(manifest.read_text())
        for output in payload.get("outputs", []):
            path = Path(output["path"])
            if not path.is_absolute():
                path = REPOSITORY / path
            if not path.is_file():
                fail(f"release output missing: {manifest}: {path}", failures)
            elif file_sha256(path) != output["sha256"]:
                fail(f"release checksum mismatch: {manifest}: {path}", failures)

    if failures:
        print(f"\nAudit failed with {len(failures)} issue(s).")
        return 1
    print(
        "Audit passed: experiment modules, driver CLIs, manuscript assets, module boundaries, "
        "cluster scripts, compatibility links, and release checksums are valid."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
