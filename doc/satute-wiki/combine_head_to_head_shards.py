#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def aggregate(detail_path, summary_path):
    groups = defaultdict(lambda: {"evaluated": 0, "informative": 0, "missing": 0})
    with open(detail_path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = (
                row["tree_case"],
                row["simulation_model"],
                row["evaluation_model"],
                row["nsites"],
                row["branch_length"],
                row["scenario"],
                row["formula"],
            )
            if row["target_found"] == "1":
                groups[key]["evaluated"] += 1
                groups[key]["informative"] += 1 if row["decision"] == "informative" else 0
            else:
                groups[key]["missing"] += 1

    with open(summary_path, "w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "tree_case",
            "simulation_model",
            "evaluation_model",
            "nsites",
            "branch_length",
            "scenario",
            "formula",
            "evaluated",
            "informative",
            "missing_split",
            "fraction_informative",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for key in sorted(groups, key=lambda x: (x[0], x[1], x[2], int(x[3]), float(x[4]), x[5], x[6])):
            group = groups[key]
            evaluated = group["evaluated"]
            fraction = group["informative"] / evaluated if evaluated else ""
            writer.writerow(
                {
                    "tree_case": key[0],
                    "simulation_model": key[1],
                    "evaluation_model": key[2],
                    "nsites": key[3],
                    "branch_length": key[4],
                    "scenario": key[5],
                    "formula": key[6],
                    "evaluated": evaluated,
                    "informative": group["informative"],
                    "missing_split": group["missing"],
                    "fraction_informative": fraction,
                }
            )


def main():
    parser = argparse.ArgumentParser(description="Combine head-to-head SatuTe shard detail TSVs.")
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    detail_path = outdir / "head_to_head_detail.tsv"
    summary_path = outdir / "head_to_head_summary.tsv"

    wrote_header = False
    with open(detail_path, "w", encoding="utf-8", newline="") as out_handle:
        writer = None
        for input_path in args.inputs:
            with open(input_path, "r", encoding="utf-8") as in_handle:
                reader = csv.DictReader(in_handle, delimiter="\t")
                if writer is None:
                    writer = csv.DictWriter(out_handle, fieldnames=reader.fieldnames, delimiter="\t")
                if not wrote_header:
                    writer.writeheader()
                    wrote_header = True
                for row in reader:
                    writer.writerow(row)

    aggregate(detail_path, summary_path)
    print(f"Detail:  {detail_path}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
