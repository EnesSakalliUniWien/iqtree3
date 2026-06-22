#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict


def load_rows(path):
    with open(path, "r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser(description="Verify paired SatuTe formula calculation invariants.")
    parser.add_argument("--detail", required=True)
    parser.add_argument("--z-tolerance", type=float, default=1e-6)
    args = parser.parse_args()

    rows = load_rows(args.detail)
    dominant_diffs = []
    for row in rows:
        if row["formula"] != "dominant" or row["target_found"] != "1" or not row["iqtree_satZ"]:
            continue
        dominant_diffs.append(abs(float(row["satZ"]) - float(row["iqtree_satZ"])))

    max_dominant_diff = max(dominant_diffs) if dominant_diffs else None
    if max_dominant_diff is None:
        raise SystemExit("No dominant rows with IQ-TREE SatuTe output were found.")
    if max_dominant_diff > args.z_tolerance:
        raise SystemExit(
            f"Dominant formula does not match IQ-TREE SatuTe: max |Zdiff|={max_dominant_diff}"
        )

    jc_by_replicate = defaultdict(dict)
    for row in rows:
        if row["model"] != "JC" or row["target_found"] != "1":
            continue
        key = (
            row["tree_case"],
            row["nsites"],
            row["branch_length"],
            row["replicate"],
            row["scenario"],
        )
        jc_by_replicate[key][row["formula"]] = float(row["satZ"])

    mismatches = []
    required = {"dominant", "all", "eigenvalue_weighted"}
    for key, values in jc_by_replicate.items():
        if not required <= values.keys():
            continue
        spread = max(values.values()) - min(values.values())
        if spread > args.z_tolerance:
            mismatches.append((key, spread))

    if mismatches:
        first_key, first_spread = mismatches[0]
        raise SystemExit(
            f"JC formulas do not collapse for {len(mismatches)} replicate/scenario rows; "
            f"first={first_key}, Z spread={first_spread}"
        )

    print(f"Rows checked: {len(rows)}")
    print(f"Dominant vs IQ-TREE max |Zdiff|: {max_dominant_diff:.6g}")
    print(f"JC collapse mismatches: {len(mismatches)}")


if __name__ == "__main__":
    main()
