#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict


def load_rows(path):
    with open(path, "r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser(description="Verify paired SatuTe calculation invariants.")
    parser.add_argument("--detail", required=True)
    parser.add_argument("--z-tolerance", type=float, default=1e-6)
    args = parser.parse_args()

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
