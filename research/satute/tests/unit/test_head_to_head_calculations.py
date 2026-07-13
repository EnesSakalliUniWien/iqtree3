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
    native_diffs = defaultdict(list)
    for row in rows:
        if row["target_found"] != "1" or not row["iqtree_satZ"]:
            continue
        native_diffs[row["formula"]].append(abs(float(row["satZ"]) - float(row["iqtree_satZ"])))

    required_native = {"dominant", "eigenvalue_weighted"}
    missing_native = sorted(required_native - set(native_diffs))
    if missing_native:
        raise SystemExit(f"No IQ-TREE SatuTe comparison rows were found for: {missing_native}")
    max_native_diffs = {formula: max(native_diffs[formula]) for formula in sorted(required_native)}
    for formula, difference in max_native_diffs.items():
        if difference > args.z_tolerance:
            raise SystemExit(
                f"{formula} does not match IQ-TREE SatuTe: max |Zdiff|={difference}"
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
    required = {"dominant", "eigenvalue_weighted"}
    for key, values in jc_by_replicate.items():
        if not required <= values.keys():
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
    for formula, difference in max_native_diffs.items():
        print(f"{formula} vs IQ-TREE max |Zdiff|: {difference:.6g}")
    print(f"JC collapse mismatches: {len(mismatches)}")


if __name__ == "__main__":
    main()
