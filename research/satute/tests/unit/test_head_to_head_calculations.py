#!/usr/bin/env python3

import argparse
import csv
import importlib.util
from collections import defaultdict
from pathlib import Path


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


def load_rows(path):
    with open(path, "r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser(description="Verify paired SatuTe calculation invariants.")
    parser.add_argument("--detail")
    parser.add_argument("--z-tolerance", type=float, default=1e-6)
    args = parser.parse_args()

    check_protein_model_contract()
    if not args.detail:
        print("Protein model contract: passed")
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
