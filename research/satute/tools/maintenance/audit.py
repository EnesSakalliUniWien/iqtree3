#!/usr/bin/env python3
"""Audit project boundaries, registry paths, and manuscript dependencies."""

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


def fail(message, failures):
    failures.append(message)
    print(f"FAIL: {message}")


def main():
    failures = []
    registry = json.loads((PROJECT / "experiments/registry.yaml").read_text())
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(PROJECT / "src") + os.pathsep + environment.get(
        "PYTHONPATH", ""
    )
    for experiment in registry["experiments"]:
        for key in ("driver", "config"):
            path = PROJECT / experiment[key]
            if not path.exists():
                fail(f"registry path missing: {experiment['id']} {key}={path}", failures)
        driver = PROJECT / experiment["driver"]
        if driver.exists():
            result = subprocess.run(
                [sys.executable, str(driver), "--help"],
                cwd=PROJECT,
                env=environment,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                text=True,
            )
            if result.returncode:
                fail(
                    f"driver CLI failed: {experiment['id']}: {result.stderr[-500:]}",
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
        "Audit passed: registry CLIs, manuscript assets, module boundaries, "
        "cluster scripts, compatibility links, and release checksums are valid."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
