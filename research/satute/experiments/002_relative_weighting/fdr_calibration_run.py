#!/usr/bin/env python3
"""Truth-labelled exact-null and mixed-null SatuTe branch-family calibration."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import math
import os
import random
import sys
import time
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]
RUN_PATH = Path(__file__).with_name("run.py")
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from satute_analysis.benchmark_contracts import FORMULAS, SCHEMA_VERSION  # noqa: E402
from satute_analysis.fdr_calibration import (  # noqa: E402
    CALIBRATION_SUMMARY_FIELDS,
    REPLICATE_FDP_FIELDS,
    TRUTH_BRANCH_FIELDS,
    aggregate_replicate_fdp,
    truth_labelled_family,
)
from satute_analysis.native_results import branch_audit_detail, read_native_rows  # noqa: E402


TIMING_FIELDS = (
    "schema_version",
    "design",
    "tree_case",
    "simulation_model",
    "nsites",
    "branch_length",
    "replicate",
    "seed",
    "truth",
    "simulation_seconds",
    "analysis_seconds",
    "total_seconds",
)


def load_run_module():
    spec = importlib.util.spec_from_file_location("relative_weighting_run", RUN_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.np = np
    return module


def calibration_seed(design, tree_case, model, nsites, branch_length, replicate, truth):
    payload = (
        f"fdr-v1:{design}:{tree_case}:{model}:{nsites}:"
        f"{float(branch_length):.12g}:{replicate}:{truth}"
    ).encode("utf-8")
    value = int.from_bytes(hashlib.blake2s(payload, digest_size=8).digest(), "big")
    return 1 + value % 2147483646


def sample_states(rng, probabilities, count):
    return rng.choice(len(probabilities), size=count, p=probabilities)


def sample_descendants(rng, parent_states, transition):
    children = np.empty_like(parent_states)
    for state in range(transition.shape[0]):
        indices = np.flatnonzero(parent_states == state)
        if indices.size:
            probabilities = np.clip(transition[state], 0.0, 1.0)
            probabilities = probabilities / probabilities.sum()
            children[indices] = rng.choice(
                transition.shape[0], size=indices.size, p=probabilities
            )
    return children


def simulate_component(start, blocked, root_states, transitions, rng):
    states = {start.id: root_states}
    sequences = {}
    stack = [(start, blocked)]
    while stack:
        node, parent = stack.pop()
        node_states = states[node.id]
        if node.name:
            sequences[node.name] = node_states
        for child, length in node.neighbors:
            if child is parent:
                continue
            transition = transitions[round(float(length), 12)]
            states[child.id] = sample_descendants(rng, node_states, transition)
            stack.append((child, node))
    return sequences


def simulate_truth_alignment(module, tree_file, model, nsites, seed, truth, output):
    base_model, gamma_shape = module.split_fixed_gamma(model)
    if gamma_shape is not None or "+I{" in base_model or "+R" in base_model:
        raise ValueError(
            "Truth-labelled calibration currently requires a homogeneous fixed DNA model; "
            "rate mixtures need shared category draws across independent components."
        )
    q, pi = module.build_q(model)
    if len(pi) != 4:
        raise ValueError("Truth-labelled calibration currently supports DNA models only")
    evals, eigenvectors, inverse = module.reversible_eigendecomposition(q, pi)
    root = module.parse_newick(tree_file)
    target_taxa = module.target_taxa_for_case(tree_file.stem)
    target_edge = module.find_edge_for_split(root, target_taxa)
    if target_edge is None:
        raise ValueError(f"Could not find target split {target_taxa} in {tree_file}")
    left, right, target_length, _side, _other_side = target_edge
    lengths = {
        round(float(length), 12)
        for node in module.all_nodes(root)
        for _other, length in node.neighbors
    }
    transitions = {
        length: module.transition_matrix(evals, eigenvectors, inverse, length)
        for length in lengths
    }
    rng = np.random.default_rng(seed)
    left_root = sample_states(rng, pi, nsites)
    if truth == "null":
        right_root = sample_states(rng, pi, nsites)
    elif truth == "alternative":
        right_root = sample_descendants(
            rng, left_root, transitions[round(float(target_length), 12)]
        )
    else:
        raise ValueError(f"Unknown truth state: {truth}")
    sequences = simulate_component(left, right, left_root, transitions, rng)
    sequences.update(simulate_component(right, left, right_root, transitions, rng))
    expected_taxa = {node.name for node in module.all_nodes(root) if node.name}
    if set(sequences) != expected_taxa:
        raise ValueError(
            f"Simulation did not cover every taxon: missing={sorted(expected_taxa - set(sequences))}"
        )
    alphabet = np.asarray(list("ACGT"))
    with output.open("w", encoding="utf-8") as handle:
        for taxon in sorted(sequences):
            handle.write(f">{taxon}\n{''.join(alphabet[sequences[taxon]])}\n")


def write_rows(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--iqtree", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--design", choices=("exact_null", "mixed_null"), required=True)
    parser.add_argument("--reps", type=int, required=True)
    parser.add_argument("--training-reps", type=int, default=0)
    parser.add_argument("--models", default="JC,GTR_PF06346")
    parser.add_argument("--site-lengths", default="100,1000")
    parser.add_argument("--branch-lengths", default="4,8")
    parser.add_argument("--tree-cases", default="five_external,sixteen_internal")
    parser.add_argument("--evonaps-branch-lengths", default="")
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    args = parser.parse_args()

    if args.reps < 1 or args.shard_count < 1:
        raise SystemExit("--reps and --shard-count must be positive")
    if not 0 <= args.shard_index < args.shard_count:
        raise SystemExit("--shard-index must satisfy 0 <= index < count")
    if args.design == "exact_null" and not 0 < args.training_reps < args.reps:
        raise SystemExit(
            "Exact-null calibration requires 0 < --training-reps < --reps"
        )
    if args.design == "mixed_null" and args.training_reps != 0:
        raise SystemExit("Mixed-null calibration requires --training-reps=0")
    if args.design == "mixed_null" and args.reps % 2:
        raise SystemExit("Mixed-null calibration requires an even --reps for a 50:50 design")
    iqtree = Path(args.iqtree)
    if not iqtree.is_file() or not os.access(iqtree, os.X_OK):
        raise SystemExit(f"Cannot execute IQ-TREE binary: {iqtree}")
    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit(f"Refusing to overwrite non-empty output directory: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    runs = outdir / "runs"
    runs.mkdir()

    module = load_run_module()
    pools = module.load_branch_length_pools(args.evonaps_branch_lengths)
    models = [
        (alias, module.resolve_model_alias(alias))
        for alias in module.parse_csv_text(args.models)
    ]
    for alias, resolved_model in models:
        if not module.is_fixed_model_spec(resolved_model):
            raise SystemExit(f"Calibration model must be fixed: {alias}={resolved_model}")
    sites = module.parse_csv_numbers(args.site_lengths, int)
    branches = module.parse_csv_numbers(args.branch_lengths, float)
    tree_cases = module.parse_csv_text(args.tree_cases)

    truth_rows = []
    replicate_rows = []
    timing_rows = []
    task_index = 0
    selected_tasks = 0
    for tree_case in tree_cases:
        for model_alias, resolved_model in models:
            for nsites in sites:
                for branch_length in branches:
                    for replicate in range(1, args.reps + 1):
                        current = task_index
                        task_index += 1
                        if current % args.shard_count != args.shard_index:
                            continue
                        selected_tasks += 1
                        if args.design == "exact_null":
                            truth = "null"
                            sample = "training" if replicate <= args.training_reps else "holdout"
                        else:
                            truth = "null" if replicate % 2 else "alternative"
                            sample = "operating"
                        seed = calibration_seed(
                            args.design,
                            tree_case,
                            model_alias,
                            nsites,
                            branch_length,
                            replicate,
                            truth,
                        )
                        case = (
                            runs
                            / tree_case
                            / f"model_{module.sanitize_model(model_alias)}"
                            / f"n{nsites}"
                            / f"b{branch_length:.2f}"
                            / f"r{replicate:05d}_{truth}"
                        )
                        case.mkdir(parents=True)
                        tree_file = case / f"{tree_case}.tree"
                        module.write_tree(
                            tree_case,
                            branch_length,
                            tree_file,
                            random.Random(seed),
                            pools,
                        )
                        started = time.perf_counter()
                        alignment = case / "alignment.fa"
                        simulation_started = time.perf_counter()
                        simulate_truth_alignment(
                            module,
                            tree_file,
                            resolved_model,
                            nsites,
                            seed,
                            truth,
                            alignment,
                        )
                        simulation_seconds = time.perf_counter() - simulation_started
                        analysis_started = time.perf_counter()
                        _tree, stat = module.run_satute(
                            str(iqtree),
                            alignment,
                            case / "true_fixed",
                            resolved_model,
                            tree_file,
                            True,
                            seed,
                        )
                        analysis_seconds = time.perf_counter() - analysis_started
                        native_rows = read_native_rows(stat)
                        branch_rows = []
                        target_taxa = module.target_taxa_for_case(tree_case)
                        for formula in FORMULAS:
                            for row in branch_audit_detail(
                                [native for native in native_rows if native["Formula"] == formula],
                                target_taxa,
                            ):
                                branch_rows.append({"formula": formula, **row})
                        base = {
                            "design": args.design,
                            "sample": sample,
                            "tree_case": tree_case,
                            "simulation_model": model_alias,
                            "evaluation_model": model_alias,
                            "nsites": nsites,
                            "branch_length": branch_length,
                            "replicate": replicate,
                            "seed": seed,
                        }
                        task_truth, task_replicate = truth_labelled_family(
                            branch_rows,
                            [target_taxa] if truth == "null" else [],
                            base,
                        )
                        truth_rows.extend(task_truth)
                        replicate_rows.extend(task_replicate)
                        timing_rows.append(
                            {
                                "schema_version": SCHEMA_VERSION,
                                "design": args.design,
                                "tree_case": tree_case,
                                "simulation_model": model_alias,
                                "nsites": nsites,
                                "branch_length": branch_length,
                                "replicate": replicate,
                                "seed": seed,
                                "truth": truth,
                                "simulation_seconds": f"{simulation_seconds:.6f}",
                                "analysis_seconds": f"{analysis_seconds:.6f}",
                                "total_seconds": f"{time.perf_counter() - started:.6f}",
                            }
                        )

    write_rows(outdir / "truth_branch_audit.tsv", TRUTH_BRANCH_FIELDS, truth_rows)
    write_rows(outdir / "replicate_fdp.tsv", REPLICATE_FDP_FIELDS, replicate_rows)
    write_rows(
        outdir / "calibration_summary.tsv",
        CALIBRATION_SUMMARY_FIELDS,
        aggregate_replicate_fdp(replicate_rows),
    )
    write_rows(outdir / "timing.tsv", TIMING_FIELDS, timing_rows)
    print(f"Shard: {args.shard_index}/{args.shard_count}")
    print(f"Tasks: {selected_tasks}/{task_index}")
    print(f"Truth rows: {len(truth_rows)}")
    print(f"Replicate rows: {len(replicate_rows)}")


if __name__ == "__main__":
    main()
