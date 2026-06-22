#!/usr/bin/env python3

import argparse
import csv
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path


PAPER_BRANCH_LENGTHS = "0.1,0.2,0.3,0.4,0.5,0.8,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,7.5,10.0"


def run(cmd):
    subprocess.run(cmd, check=True)


def parse_csv_numbers(text, cast=float):
    return [cast(value) for value in text.split(",") if value.strip()]


def pair(prefix, i, j, length):
    return f"({prefix}{i}:{length:.6g},{prefix}{j}:{length:.6g}):{length:.6g}"


def eight_taxon_subtree(prefix, length):
    left = f"({pair(prefix, 1, 2, length)},{pair(prefix, 3, 4, length)}):{length:.6g}"
    right = f"({pair(prefix, 5, 6, length)},{pair(prefix, 7, 8, length)}):{length:.6g}"
    return f"({left},{right})"


def five_taxon_tree(branch_length):
    return f"(A:{branch_length:.6g},(B:0.2,C:0.2):0.2,(D:0.2,E:0.2):0.2);\n"


def sixteen_taxon_tree(branch_length, subtree_length):
    a_side = eight_taxon_subtree("A", subtree_length)
    b_side = eight_taxon_subtree("B", subtree_length)
    return f"({a_side}:{branch_length:.6g},{b_side}:0.0);\n"


def parse_sat_stat(path):
    header = None
    rows = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("#") or not line:
                continue
            fields = line.split("\t")
            if fields[0] == "ID":
                header = fields
                continue
            if header and fields[0].isdigit():
                rows.append(dict(zip(header, fields)))
    return rows


def target_row(stat_file, split):
    for row in parse_sat_stat(stat_file):
        if row.get("Split") == split:
            return row
    return None


def write_row(writer, tree_case, nsites, branch_length, rep, scenario, row, decision_column):
    if row is None:
        writer.writerow(
            {
                "tree_case": tree_case,
                "nsites": nsites,
                "branch_length": branch_length,
                "replicate": rep,
                "scenario": scenario,
                "target_found": 0,
                "tested_length": "",
                "satZ": "",
                "satP": "",
                "decision": "",
                "decision_bonferroni": "",
                "informative": "",
            }
        )
        return

    decision = row[decision_column]
    writer.writerow(
        {
            "tree_case": tree_case,
            "nsites": nsites,
            "branch_length": branch_length,
            "replicate": rep,
            "scenario": scenario,
            "target_found": 1,
            "tested_length": row["Length"],
            "satZ": row["satZ"],
            "satP": row["satP"],
            "decision": row["Decision"],
            "decision_bonferroni": row["DecisionTaxonBonf"],
            "informative": 1 if decision == "informative" else 0,
        }
    )


def run_satute(iqtree, alignment, prefix, split, tree=None, fixed_lengths=False):
    cmd = [
        iqtree,
        "-s",
        str(alignment),
    ]
    if tree is not None:
        cmd.extend(["-te", str(tree)])
    cmd.extend(
        [
            "-m",
            "JC",
            "--satute",
            "--prefix",
            str(prefix),
            "-T",
            "1",
            "--redo",
            "--quiet",
        ]
    )
    if fixed_lengths:
        cmd.append("-blfix")
    run(cmd)
    return target_row(str(prefix) + ".sat.stat", split)


def run_case(iqtree, outdir, writer, tree_case, nsites, branch_length, rep, seed):
    case_dir = outdir / "runs" / tree_case / f"n{nsites}" / f"b{branch_length:.2f}" / f"r{rep:03d}"
    case_dir.mkdir(parents=True, exist_ok=True)
    tree_file = case_dir / "true.tree"

    if tree_case == "five_external":
        split = "A"
        tree_file.write_text(five_taxon_tree(branch_length), encoding="utf-8")
    elif tree_case == "sixteen_internal":
        split = ",".join([f"A{i}" for i in range(1, 9)])
        tree_file.write_text(sixteen_taxon_tree(branch_length, 0.1), encoding="utf-8")
    else:
        raise ValueError(tree_case)

    sim_prefix = case_dir / "sim"
    run(
        [
            iqtree,
            "--alisim",
            str(sim_prefix),
            "-t",
            str(tree_file),
            "-m",
            "JC",
            "--length",
            str(nsites),
            "--seed",
            str(seed),
            "-af",
            "fasta",
            "--quiet",
        ]
    )
    alignment = Path(str(sim_prefix) + ".fa")

    row = run_satute(iqtree, alignment, case_dir / "true_fixed", split, tree_file, fixed_lengths=True)
    write_row(writer, tree_case, nsites, branch_length, rep, "true_tree_fixed_lengths", row, "Decision")

    row = run_satute(iqtree, alignment, case_dir / "true_topology_ml_lengths", split, tree_file, fixed_lengths=False)
    write_row(writer, tree_case, nsites, branch_length, rep, "true_topology_ml_lengths", row, "Decision")

    row = run_satute(iqtree, alignment, case_dir / "ml_tree", split, tree=None, fixed_lengths=False)
    write_row(writer, tree_case, nsites, branch_length, rep, "ml_tree_unadjusted", row, "Decision")
    write_row(writer, tree_case, nsites, branch_length, rep, "ml_tree_bonferroni", row, "DecisionTaxonBonf")


def aggregate(detail_path, summary_path):
    groups = defaultdict(lambda: {"evaluated": 0, "informative": 0, "missing": 0})
    with open(detail_path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = (row["tree_case"], row["nsites"], row["branch_length"], row["scenario"])
            if row["target_found"] == "1":
                groups[key]["evaluated"] += 1
                groups[key]["informative"] += int(row["informative"])
            else:
                groups[key]["missing"] += 1

    with open(summary_path, "w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "tree_case",
            "nsites",
            "branch_length",
            "scenario",
            "evaluated",
            "informative",
            "missing_split",
            "fraction_informative",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for key in sorted(groups, key=lambda item: (item[0], int(item[1]), float(item[2]), item[3])):
            group = groups[key]
            fraction = group["informative"] / group["evaluated"] if group["evaluated"] else ""
            writer.writerow(
                {
                    "tree_case": key[0],
                    "nsites": key[1],
                    "branch_length": key[2],
                    "scenario": key[3],
                    "evaluated": group["evaluated"],
                    "informative": group["informative"],
                    "missing_split": group["missing"],
                    "fraction_informative": fraction,
                }
            )


def main():
    parser = argparse.ArgumentParser(
        description="Reduced or full reconstruction of the SatuTe Fig. 2 JC simulation design."
    )
    parser.add_argument("--iqtree", default="./build/iqtree3")
    parser.add_argument("--outdir", default="/tmp/iqtree-satute-paper-fig2-reconstruction")
    parser.add_argument("--reps", type=int, default=3)
    parser.add_argument("--site-lengths", default="100,1000")
    parser.add_argument("--branch-lengths", default="0.1,1.0,5.0")
    parser.add_argument("--paper-grid", action="store_true", help="Use the paper branch-length grid and site lengths 100,1000,10000.")
    parser.add_argument("--tree-cases", default="five_external,sixteen_internal")
    args = parser.parse_args()

    iqtree = str(Path(args.iqtree))
    if not Path(iqtree).exists():
        raise SystemExit(f"IQ-TREE binary not found: {iqtree}")

    outdir = Path(args.outdir)
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True)

    branch_lengths = parse_csv_numbers(PAPER_BRANCH_LENGTHS if args.paper_grid else args.branch_lengths, float)
    site_lengths = parse_csv_numbers("100,1000,10000" if args.paper_grid else args.site_lengths, int)
    tree_cases = [value.strip() for value in args.tree_cases.split(",") if value.strip()]

    detail_path = outdir / "fig2_reconstruction_detail.tsv"
    summary_path = outdir / "fig2_reconstruction_summary.tsv"
    with open(detail_path, "w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "tree_case",
            "nsites",
            "branch_length",
            "replicate",
            "scenario",
            "target_found",
            "tested_length",
            "satZ",
            "satP",
            "decision",
            "decision_bonferroni",
            "informative",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for tree_case in tree_cases:
            for nsites in site_lengths:
                for branch_length in branch_lengths:
                    for rep in range(1, args.reps + 1):
                        seed = 100000 + rep + nsites * 10 + int(branch_length * 1000)
                        run_case(iqtree, outdir, writer, tree_case, nsites, branch_length, rep, seed)

    aggregate(detail_path, summary_path)
    print(f"Detail:  {detail_path}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
