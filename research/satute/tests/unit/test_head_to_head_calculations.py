#!/usr/bin/env python3

import argparse
import csv
import importlib.util
import random
import tempfile
from collections import defaultdict
from pathlib import Path
from unittest import mock


RUN_PATH = Path(__file__).resolve().parents[2] / "experiments" / "002_relative_weighting" / "run.py"


def load_run_module():
    spec = importlib.util.spec_from_file_location("relative_weighting_run", RUN_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def check_protein_model_contract():
    module = load_run_module()
    expected = {
        "LG_G4": "LG+G4{0.5}",
        "WAG_G4": "WAG+G4{0.5}",
        "JTT_G4": "JTT+G4{0.5}",
        "Q.PFAM_G4": "Q.pfam+G4{0.5}",
    }
    observed = {alias: module.resolve_model_alias(alias) for alias in expected}
    if observed != expected:
        raise SystemExit(f"Protein alias contract mismatch: {observed}")
    command = module.seqgen_command("seq-gen", "WAG+G4{0.5}", 100, 7)
    if command != ["seq-gen", "-mWAG", "-a", "0.5", "-g", "4", "-l", "100", "-n", "1", "-z", "7"]:
        raise SystemExit(f"Unexpected Seq-Gen WAG+G4 command: {command}")
    if not module.is_fixed_model_spec("Q.PFAM_G4"):
        raise SystemExit("Q.PFAM_G4 must be treated as a fixed model")
    gtr20 = "GTR20{" + ",".join(["1"] * 189) + "}+F{" + ",".join(["0.05"] * 20) + "}+G4{0.5}"
    rates, frequencies = module.parse_model(gtr20)
    if len(rates) != 190 or rates[-1] != 1.0 or len(frequencies) != 20:
        raise SystemExit("IQ-TREE's 189-parameter GTR20 syntax was not expanded for Seq-Gen")


def check_fixed_dna_mixture_contract():
    module = load_run_module()
    base = module.GTR_PF06346_MODEL
    expected = {
        "GTR_PF06346_G4": base + "+G4{0.5}",
        "GTR_PF06346_I_G4": base + "+I{0.1}+G4{0.5}",
    }
    observed = {alias: module.resolve_model_alias(alias) for alias in expected}
    if observed != expected:
        raise SystemExit(f"Fixed DNA mixture alias contract mismatch: {observed}")
    for alias in expected:
        if not module.is_fixed_model_spec(alias):
            raise SystemExit(f"{alias} must be treated as a fixed evaluation model")
    gamma_command = module.seqgen_command("seq-gen", expected["GTR_PF06346_G4"], 100, 7)
    if gamma_command[-10:] != [
        "-a", "0.5", "-g", "4", "-l", "100", "-n", "1", "-z", "7"
    ]:
        raise SystemExit(f"Unexpected Seq-Gen fixed GTR+G4 command: {gamma_command}")
    try:
        module.seqgen_command("seq-gen", expected["GTR_PF06346_I_G4"], 100, 7)
    except ValueError as exc:
        if "--simulator alisim" not in str(exc):
            raise SystemExit(f"Unexpected +I backend error: {exc}") from exc
    else:
        raise SystemExit("Seq-Gen must reject +I models rather than silently dropping invariant sites")
    pairs = module.resolve_requested_model_pairs(
        "GTR_PF06346_G4:GTR_PF06346_G4,GTR_PF06346_I_G4:GTR_PF06346_I_G4",
        "",
        "",
    )
    description = module.format_resolved_model_pairs(pairs)
    if expected["GTR_PF06346_G4"] not in description or expected["GTR_PF06346_I_G4"] not in description:
        raise SystemExit(f"Resolved model-pair description lost fixed parameters: {description}")


def check_unrooted_tree_contract():
    module = load_run_module()
    with tempfile.TemporaryDirectory() as temporary:
        temporary = Path(temporary)
        cases = (
            (
                "five_external",
                module.five_taxon_tree(0.7),
            ),
            (
                "sixteen_internal",
                module.sixteen_taxon_tree(0.7, random.Random(17), {"internal": [], "external": []}),
            ),
        )
        for tree_case, newick in cases:
            path = temporary / f"{tree_case}.tree"
            path.write_text(newick, encoding="utf-8")
            root = module.parse_newick(path)
            if len(root.neighbors) != 3:
                raise SystemExit(f"{tree_case} benchmark tree is not explicitly unrooted: {newick.strip()}")
            edge = module.find_edge_for_split(root, module.target_taxa_for_case(tree_case))
            if edge is None or abs(edge[2] - 0.7) > 1e-12:
                raise SystemExit(f"{tree_case} focal branch was not preserved: {edge}")


def check_simulation_cache_contract():
    module = load_run_module()
    cache = {}
    with tempfile.TemporaryDirectory() as temporary:
        temporary = Path(temporary)

        def fake_simulator(_iqtree, _seqgen, _indelible, _backend, _tree, _model, _nsites, _seed, prefix):
            alignment = Path(str(prefix) + ".fa")
            alignment.write_text(">A1\nAAAA\n>A2\nAAAA\n>A3\nAAAA\n>A4\nAAAA\n>B\nAAAA\n", encoding="utf-8")
            return alignment

        with mock.patch.object(module, "simulate_alignment", side_effect=fake_simulator) as simulator:
            positional = (
                "iqtree", "seq-gen", "indelible", "alisim", {"internal": [], "external": []},
                temporary, cache,
            )
            # Exercise only the cache-creation logic through run_case by
            # replacing the expensive native-analysis stages with no-ops.
            with mock.patch.object(module, "selected_scenarios", return_value=[]):
                writer = mock.Mock()
                timing_writer = mock.Mock()
                trailing = (
                    writer, writer, writer, timing_writer, "five_external", module.GTR_PF06346_MODEL,
                    module.GTR_PF06346_MODEL, 4, 0.5, 1, 17, 0.05, "fig2",
                )
                module.run_case(*positional, *trailing)
                module.run_case(*positional, *trailing)
        if simulator.call_count != 1:
            raise SystemExit(f"Simulation cache did not reuse the alignment: calls={simulator.call_count}")
        timing_rows = [call.args[0] for call in timing_writer.writerow.call_args_list]
        if [row["simulation_cache_hit"] for row in timing_rows] != [0, 1]:
            raise SystemExit(f"Unexpected cache-hit timing flags: {timing_rows}")


def load_rows(path):
    with open(path, "r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser(description="Verify paired SatuTe calculation invariants.")
    parser.add_argument("--detail")
    parser.add_argument("--z-tolerance", type=float, default=1e-6)
    args = parser.parse_args()

    check_protein_model_contract()
    check_fixed_dna_mixture_contract()
    check_unrooted_tree_contract()
    check_simulation_cache_contract()
    if not args.detail:
        print("Protein and unrooted-tree contracts: passed")
        return

    rows = load_rows(args.detail)
    required_formulas = {"dominant", "eigenvalue_weighted"}
    unexpected_formulas = sorted({row["formula"] for row in rows} - required_formulas)
    if unexpected_formulas:
        raise SystemExit(f"Unexpected benchmark formulas: {unexpected_formulas}")
    non_native = [row for row in rows if row.get("implementation") != "iqtree_native"]
    if non_native:
        raise SystemExit(f"Found {len(non_native)} rows not sourced from native IQ-TREE output")

    formulas_by_case = defaultdict(set)
    for row in rows:
        key = (
            row["tree_case"],
            row["simulation_model"],
            row["evaluation_model"],
            row["nsites"],
            row["branch_length"],
            row["replicate"],
            row["scenario"],
        )
        formulas_by_case[key].add(row["formula"])
    incomplete = [
        (key, formulas)
        for key, formulas in formulas_by_case.items()
        if formulas != required_formulas
    ]
    if incomplete:
        key, formulas = incomplete[0]
        raise SystemExit(
            f"Native formula pair is incomplete for {len(incomplete)} case/scenario groups; "
            f"first={key}, formulas={sorted(formulas)}"
        )

    jc_by_replicate = defaultdict(dict)
    for row in rows:
        if row["evaluation_model"] != "JC" or row["target_found"] != "1":
            continue
        key = (
            row["tree_case"],
            row["simulation_model"],
            row["nsites"],
            row["branch_length"],
            row["replicate"],
            row["scenario"],
        )
        jc_by_replicate[key][row["formula"]] = float(row["satZ"])

    mismatches = []
    for key, values in jc_by_replicate.items():
        if not required_formulas <= values.keys():
            continue
        spread = abs(values["dominant"] - values["eigenvalue_weighted"])
        if spread > args.z_tolerance:
            mismatches.append((key, spread))

    if mismatches:
        first_key, first_spread = mismatches[0]
        raise SystemExit(
            f"JC dominant and eigenvalue-weighted statistics do not collapse for {len(mismatches)} replicate/scenario rows; "
            f"first={first_key}, Z spread={first_spread}"
        )

    print(f"Rows checked: {len(rows)}")
    print(f"Native formula pairs checked: {len(formulas_by_case)}")
    print(f"JC collapse mismatches: {len(mismatches)}")


if __name__ == "__main__":
    main()
