"""Single public command interface for the SatuTe research workspace."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

from .cluster_campaign import CampaignError, load_campaign, print_plan, submit_campaign
from .project import PROJECT_ROOT


REPOSITORY_ROOT = PROJECT_ROOT.parents[1]
ENTRY_POINT = PROJECT_ROOT / "satute"
VENV_ROOT = PROJECT_ROOT / ".venv"
IQTREE = REPOSITORY_ROOT / "build" / "iqtree3"
DEVELOPER_ADAPTERS = {
    "biological-analysis": "experiments/005_biological_windows/analyze_dataset.py",
    "formula-comparison": "tests/native/compare_formulas.py",
    "invariant-full-figure": "experiments/003_invariant_mixture/render_full.py",
    "invariant-shift-figure": "experiments/003_invariant_mixture/render_shift.py",
    "invariant-simulation": "experiments/003_invariant_mixture/run.py",
    "paper-head-to-head": "experiments/002_relative_weighting/run.py",
    "reference": "tests/native/satute_reference.py",
    "reproduce-scientific-results": "tools/maintenance/reproduce_scientific_results.py",
    "saturation-null-check": "tests/exact/test_saturation_null.py",
    "saturation-simulation": "experiments/004_saturation_diagnostics/run.py",
    "validate-weighted-formula": "tests/native/validate_weighted_formula.py",
}


def run_python(relative_script: str, *arguments: object) -> None:
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(PROJECT_ROOT / "src")
    subprocess.run(
        [
            sys.executable,
            str(PROJECT_ROOT / relative_script),
            *(str(argument) for argument in arguments),
        ],
        cwd=PROJECT_ROOT,
        env=environment,
        check=True,
    )


def discover_virtual_environments() -> list[str]:
    names = {".venv", "venv", "env", "virtualenv"}
    environments: list[str] = []
    for directory in REPOSITORY_ROOT.rglob("*"):
        if not directory.is_dir() or directory.name not in names:
            continue
        if any(part in {"artifacts", "build", ".git"} for part in directory.parts):
            continue
        if (directory / "pyvenv.cfg").is_file():
            environments.append(str(directory))
    return sorted(environments)


def doctor(as_json: bool) -> int:
    expected_python = VENV_ROOT / "bin" / "python"
    report = {
        "entry_point": str(ENTRY_POINT),
        "python": str(expected_python),
        "runtime_python": str(Path(sys.executable).resolve()),
        "runtime_prefix": str(Path(sys.prefix).resolve()),
        "expected_python": str(expected_python),
        "virtual_environments": discover_virtual_environments(),
        "iqtree": str(IQTREE),
        "iqtree_available": IQTREE.is_file(),
    }
    healthy = (
        Path(sys.prefix).resolve() == VENV_ROOT.resolve()
        and report["virtual_environments"] == [str(VENV_ROOT)]
    )
    if as_json:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        print(f"Entry point: {report['entry_point']}")
        print(f"Python: {report['python']}")
        print(f"Virtual environment: {VENV_ROOT}")
        print(f"IQ-TREE: {IQTREE} ({'available' if report['iqtree_available'] else 'missing'})")
        print(f"Workspace environment: {'healthy' if healthy else 'inconsistent'}")
    return 0 if healthy else 1


def data_command(args: argparse.Namespace) -> None:
    command: list[object] = []
    if args.map:
        command.append("--map")
    else:
        command.extend(("--profile", args.profile))
    for selection in args.only:
        command.extend(("--only", selection))
    if args.plan:
        command.append("--plan")
    if args.skip_rate_category_check:
        command.append("--skip-rate-category-check")
    if args.figures:
        command.extend(("--figures", "--figure-mode", args.figure_mode))
    run_python("generate_experiment_data.py", *command)


def experiment_command(args: argparse.Namespace) -> None:
    command: list[object] = (
        ["--experiment", args.experiment, "--profile", args.profile]
    )
    if args.study:
        command.extend(("--study", args.study))
    if args.figures:
        command.append("--figures")
    run_python("tools/maintenance/run_experiment.py", *command)


def publication_figures_command(args: argparse.Namespace) -> None:
    command: list[object] = [
        "--mode",
        "index-existing" if args.index_existing else "render",
    ]
    if args.plan:
        command.append("--plan")
    run_python("tools/maintenance/generate_publication_power_figures.py", *command)


def cluster_command(args: argparse.Namespace) -> int:
    try:
        campaign = load_campaign(Path(args.campaign))
        if args.cluster_action == "plan":
            print_plan(campaign, args.only)
        elif args.cluster_action == "submit":
            submit_campaign(campaign, args.only)
        else:
            raise AssertionError(args.cluster_action)
    except (CampaignError, OSError, subprocess.CalledProcessError) as error:
        print(f"Cluster campaign error: {error}", file=sys.stderr)
        return 2
    return 0


def run_check(check_name: str) -> None:
    unit_checks = (
        "tests/unit/test_experiment_catalog.py",
        "tests/unit/test_workspace_interface.py",
        "tests/unit/test_cluster_campaign.py",
        "tests/unit/test_head_to_head_calculations.py",
        "tests/unit/test_benchmark_contracts.py",
        "tests/unit/test_publication_figure_contracts.py",
        "tests/unit/test_multiple_testing.py",
        "tests/unit/test_fdr_calibration.py",
        "tests/unit/test_fdr_postprocess.py",
    )
    if check_name == "all":
        for script in unit_checks:
            run_python(script)
        run_python("tests/exact/test_weighted_likelihood.py")
        run_python("tests/native/test_satute_mixture_null.py", IQTREE, "artifacts/local/check/mixture-null")
        run_python(
            "tests/exact/test_saturation_null.py",
            "--replicates",
            2000,
            "--output",
            "artifacts/local/check/exact-null.tsv",
        )
        environment = dict(os.environ)
        environment["PYTHON_BIN"] = sys.executable
        subprocess.run(
            [
                "bash",
                "tests/native/run_satute_phase1_checks.sh",
                str(IQTREE),
                "artifacts/local/check/native",
            ],
            cwd=PROJECT_ROOT,
            env=environment,
            check=True,
        )
        run_python("tools/maintenance/audit.py")
        return
    if check_name == "power":
        for script in unit_checks[2:]:
            run_python(script)
        run_python("tests/exact/test_weighted_likelihood.py")
        run_python("tests/native/test_satute_mixture_null.py", IQTREE, "artifacts/local/check/mixture-null")
        run_python(
            "tests/native/test_satute_rate_categories.py",
            IQTREE,
            "artifacts/local/check/rate-categories",
        )
        return
    if check_name == "figures":
        for script in (
            "tools/maintenance/audit_power_plots.py",
            "tools/maintenance/verify_figure_integration.py",
            "tools/maintenance/journal_artwork_gate.py",
            "tools/maintenance/verify_equation_integration.py",
        ):
            run_python(script)
        return
    if check_name == "publication-figures":
        run_python("tests/unit/test_publication_figure_contracts.py")
        run_python("tools/maintenance/audit_power_plots.py")
        return
    if check_name == "rate-categories":
        run_python(
            "tests/native/test_satute_rate_categories.py",
            IQTREE,
            "artifacts/local/check/rate-categories",
        )
        return
    if check_name == "audit":
        run_python("tools/maintenance/audit.py")
        return
    raise AssertionError(check_name)


def manuscript_command() -> None:
    build = PROJECT_ROOT / "manuscript" / "build"
    build.mkdir(parents=True, exist_ok=True)
    for source in ("main.tex", "supplementary.tex"):
        subprocess.run(
            [
                "latexmk",
                "-norc",
                "-pdf",
                "-interaction=nonstopmode",
                "-halt-on-error",
                "-outdir=build",
                source,
            ],
            cwd=PROJECT_ROOT / "manuscript",
            check=True,
        )
    shutil.copy2(build / "main.pdf", PROJECT_ROOT / "manuscript" / "manuscript.pdf")
    shutil.copy2(
        build / "supplementary.pdf",
        PROJECT_ROOT / "manuscript" / "supplementary-information.pdf",
    )


LOCAL_ARTIFACTS = PROJECT_ROOT / "artifacts" / "local"

# Removed by a bare `clean`: each is recreated in full by one documented
# command, so deleting it costs only the time to re-run that command.
REGENERABLE_TARGETS = (
    (PROJECT_ROOT / "manuscript" / "build", "./satute manuscript"),
    (LOCAL_ARTIFACTS / "smoke", "./satute data --profile smoke"),
    (LOCAL_ARTIFACTS / "check", "./satute check all"),
    (PROJECT_ROOT / ".pytest_cache", "./satute check all"),
)


def directory_bytes(path: Path) -> int:
    if path.is_file():
        return path.stat().st_size
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def format_bytes(count: int) -> str:
    size = float(count)
    for unit in ("B", "KB", "MB", "GB"):
        if size < 1024.0 or unit == "GB":
            return f"{size:.0f} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1024.0
    return f"{size:.1f} GB"


def declared_output_directories() -> dict[Path, str]:
    """Output directories some experiment manifest claims, and who claims them."""
    declared: dict[Path, str] = {}
    for manifest in sorted(PROJECT_ROOT.glob("experiments/*/experiment.json")):
        try:
            data = json.loads(manifest.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        output = data.get("config", {}).get("output")
        if isinstance(output, str) and output.strip():
            declared[(PROJECT_ROOT / output).resolve()] = data.get("id", manifest.parent.name)
    return declared


def report_local_artifacts() -> None:
    """Inventory artifacts/local so it is clear what is safe to delete.

    Anything a manifest declares is rebuilt by `./satute experiment <id>`.
    Anything else is a one-off exploratory run that no command recreates, so
    `clean` never touches it and removing it has to be asked for by name.
    """
    if not LOCAL_ARTIFACTS.is_dir():
        print(f"{LOCAL_ARTIFACTS} does not exist")
        return
    declared = declared_output_directories()
    regenerable = {path.resolve() for path, _ in REGENERABLE_TARGETS}

    entries = sorted(LOCAL_ARTIFACTS.iterdir(), key=lambda p: -directory_bytes(p))
    total = 0
    print(f"{'entry':<44} {'size':>9}  provenance")
    for entry in entries:
        size = directory_bytes(entry)
        total += size
        resolved = entry.resolve()
        if resolved in regenerable:
            provenance = "regenerable, removed by ./satute clean"
        elif resolved in declared:
            provenance = f"declared by {declared[resolved]}"
        else:
            provenance = "ad hoc -- no command recreates this"
        print(f"{entry.name:<44} {format_bytes(size):>9}  {provenance}")
    print(f"{'total':<44} {format_bytes(total):>9}")
    print()
    print("Remove an ad hoc entry with:  ./satute clean --purge <entry> [...]")


def purge_local_artifacts(names: list[str]) -> int:
    """Delete named entries under artifacts/local. Refuses anything else."""
    failures = 0
    for name in names:
        target = (LOCAL_ARTIFACTS / name).resolve()
        try:
            target.relative_to(LOCAL_ARTIFACTS.resolve())
        except ValueError:
            print(f"refusing {name}: outside {LOCAL_ARTIFACTS}")
            failures += 1
            continue
        if target == LOCAL_ARTIFACTS.resolve():
            print("refusing to remove artifacts/local itself; name the entries")
            failures += 1
            continue
        if not target.exists():
            print(f"skipping {name}: not present")
            continue
        size = directory_bytes(target)
        if target.is_dir():
            shutil.rmtree(target)
        else:
            target.unlink()
        print(f"removed {name} ({format_bytes(size)})")
    return failures


def clean_command(dry_run: bool = False, purge: list[str] | None = None) -> int:
    if purge:
        return 1 if purge_local_artifacts(purge) else 0

    reclaimed = 0
    for target, rebuild in REGENERABLE_TARGETS:
        if not target.is_dir():
            continue
        size = directory_bytes(target)
        reclaimed += size
        label = target.relative_to(PROJECT_ROOT)
        if dry_run:
            print(f"would remove {label} ({format_bytes(size)}) -- rebuild: {rebuild}")
        else:
            shutil.rmtree(target)
            print(f"removed {label} ({format_bytes(size)})")

    caches = [
        cache for cache in sorted(PROJECT_ROOT.rglob("__pycache__"), reverse=True)
        if cache.is_dir() and ".venv" not in cache.parts
    ]
    if caches:
        cache_bytes = sum(directory_bytes(cache) for cache in caches)
        reclaimed += cache_bytes
        if dry_run:
            print(f"would remove {len(caches)} __pycache__ directories ({format_bytes(cache_bytes)})")
        else:
            for cache in caches:
                shutil.rmtree(cache)
            print(f"removed {len(caches)} __pycache__ directories ({format_bytes(cache_bytes)})")

    print(f"{'would reclaim' if dry_run else 'reclaimed'} {format_bytes(reclaimed)}")
    print("Run ./satute clean --report to inventory the rest of artifacts/local.")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="./satute",
        description="Single interface for SatuTe experiments, figures, checks, and manuscript.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("bootstrap", help="Create/update the sole .venv environment.")

    doctor_parser = subparsers.add_parser("doctor", help="Verify entry-point and environment consistency.")
    doctor_parser.add_argument("--json", action="store_true")

    subparsers.add_parser("status", help="Report experiment and artifact status.")

    data_parser = subparsers.add_parser("data", help="Plan or run managed experiment data.")
    data_parser.add_argument("--profile", choices=("smoke", "full"), default="smoke")
    data_parser.add_argument("--map", action="store_true")
    data_parser.add_argument("--only", action="append", default=[])
    data_parser.add_argument("--plan", action="store_true")
    data_parser.add_argument("--skip-rate-category-check", action="store_true")
    data_parser.add_argument("--figures", action="store_true")
    data_parser.add_argument(
        "--figure-mode", choices=("render", "index-existing"), default="render"
    )

    experiment_parser = subparsers.add_parser("experiment", help="Run one registered experiment study.")
    experiment_parser.add_argument("experiment")
    experiment_parser.add_argument("--study")
    experiment_parser.add_argument("--profile", choices=("smoke", "full"), default="smoke")
    experiment_parser.add_argument("--figures", action="store_true")

    figures_parser = subparsers.add_parser(
        "publication-figures", help="Build the canonical publication figure suite."
    )
    figures_parser.add_argument("--plan", action="store_true")
    figures_parser.add_argument("--index-existing", action="store_true")

    cluster_parser = subparsers.add_parser(
        "cluster",
        help="Plan or submit the versioned LiSC recalculation campaign.",
    )
    cluster_subparsers = cluster_parser.add_subparsers(
        dest="cluster_action", required=True
    )
    for action in ("plan", "submit"):
        action_parser = cluster_subparsers.add_parser(action)
        action_parser.add_argument(
            "--campaign",
            default=str(PROJECT_ROOT / "campaigns/2026_expanded_recalculation.json"),
        )
        action_parser.add_argument(
            "--only",
            action="append",
            default=[],
            help="Limit the operation to this campaign run ID; repeat as needed.",
        )

    check_parser = subparsers.add_parser("check", help="Run a reproducible validation suite.")
    check_parser.add_argument(
        "suite",
        nargs="?",
        choices=("all", "power", "figures", "publication-figures", "rate-categories", "audit"),
        default="all",
    )

    subparsers.add_parser("manuscript", help="Compile the main and supplementary manuscripts.")
    subparsers.add_parser("workflow-figure", help="Render the LaTeX workflow figure.")

    release_parser = subparsers.add_parser("release", help="Create a reviewed experiment release.")
    release_parser.add_argument("experiment")
    release_parser.add_argument("version")

    developer_parser = subparsers.add_parser(
        "developer",
        help="Run an allow-listed low-level adapter through the canonical environment.",
    )
    developer_parser.add_argument("adapter", choices=tuple(DEVELOPER_ADAPTERS))
    developer_parser.add_argument("arguments", nargs=argparse.REMAINDER)

    clean_parser = subparsers.add_parser(
        "clean", help="Remove reproducible local build and smoke outputs."
    )
    clean_parser.add_argument(
        "--dry-run", action="store_true", help="Show what would be removed."
    )
    clean_parser.add_argument(
        "--report",
        action="store_true",
        help="Inventory artifacts/local with sizes and whether each entry is regenerable.",
    )
    clean_parser.add_argument(
        "--purge",
        nargs="+",
        metavar="ENTRY",
        default=None,
        help="Delete named artifacts/local entries. For ad hoc runs no command recreates.",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    if args.command == "bootstrap":
        parser.error("bootstrap is handled by the ./satute launcher")
    if args.command == "doctor":
        return doctor(args.json)
    if args.command == "status":
        run_python("tools/maintenance/status.py")
    elif args.command == "data":
        data_command(args)
    elif args.command == "experiment":
        experiment_command(args)
    elif args.command == "publication-figures":
        publication_figures_command(args)
    elif args.command == "cluster":
        return cluster_command(args)
    elif args.command == "check":
        run_check(args.suite)
    elif args.command == "manuscript":
        manuscript_command()
    elif args.command == "workflow-figure":
        run_python("tools/manuscript/render_workflow_latex.py")
    elif args.command == "release":
        run_python(
            "tools/maintenance/create_release_manifest.py",
            "--experiment",
            args.experiment,
            "--version",
            args.version,
        )
    elif args.command == "developer":
        arguments = args.arguments[1:] if args.arguments[:1] == ["--"] else args.arguments
        run_python(DEVELOPER_ADAPTERS[args.adapter], *arguments)
    elif args.command == "clean":
        if args.report:
            report_local_artifacts()
            return 0
        return clean_command(dry_run=args.dry_run, purge=args.purge)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
