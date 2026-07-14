#!/usr/bin/env python3

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


MODEL = "GTR{1,2,1,1,2,1}+F{0.30,0.20,0.20,0.30}+G4{0.5}"


def run(cmd):
    subprocess.run(cmd, check=True)


def branch_id_for_split(path, split):
    header = None
    with Path(path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[0] == "ID":
                header = fields
                continue
            if header is None or not fields[0].isdigit():
                continue
            row = dict(zip(header, fields))
            if row.get("Split") == split and row.get("Formula") == "dominant" and row.get("RateCategory") == "pooled":
                return row["ID"]
    raise ValueError(f"Could not find split {split} in {path}")


def main():
    parser = argparse.ArgumentParser(
        description="Verify the generalized SatuTe Python reference on a larger fixed tree."
    )
    parser.add_argument("iqtree", help="Path to iqtree3")
    parser.add_argument("outdir", nargs="?", default="/tmp/iqtree-satute-reference-general")
    args = parser.parse_args()

    iqtree = Path(args.iqtree)
    if not iqtree.exists():
        raise SystemExit(f"Cannot find IQ-TREE binary: {iqtree}")

    outdir = Path(args.outdir)
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True)

    tree = outdir / "six_taxon.tree"
    tree.write_text(
        "((A:0.05,B:0.05):1.4,C:0.08,(D:0.05,(E:0.05,F:0.05):0.12):0.45);\n",
        encoding="utf-8",
    )

    sim_prefix = outdir / "six_taxon"
    sat_prefix = outdir / "six_taxon_sat"
    run(
        [
            str(iqtree),
            "--alisim",
            str(sim_prefix),
            "-t",
            str(tree),
            "-m",
            MODEL,
            "--length",
            "1500",
            "--seed",
            "612",
            "-af",
            "fasta",
            "--quiet",
        ]
    )
    run(
        [
            str(iqtree),
            "-s",
            str(sim_prefix) + ".fa",
            "-te",
            str(tree),
            "-m",
            MODEL,
            "--satute",
            "--rate",
            "--prefix",
            str(sat_prefix),
            "-T",
            "1",
            "--redo",
            "--quiet",
        ]
    )

    reference = Path(__file__).with_name("satute_reference.py")
    reference_out = outdir / "six_taxon_reference.tsv"
    sat_stat = str(sat_prefix) + ".sat.stat"
    run(
        [
            sys.executable,
            str(reference),
            "--alignment",
            str(sim_prefix) + ".fa",
            "--tree",
            str(sat_prefix) + ".sat.tree",
            "--model",
            MODEL,
            "--split",
            "A,B",
            "--rate-file",
            str(sat_prefix) + ".rate",
            "--iqtree-report",
            str(sat_prefix) + ".iqtree",
            "--out",
            str(reference_out),
            "--compare-sat-stat",
            sat_stat,
        ]
    )

    branch_id = branch_id_for_split(sat_stat, "A,B")
    branch_reference_out = outdir / "six_taxon_reference_branch_id.tsv"
    run(
        [
            sys.executable,
            str(reference),
            "--alignment",
            str(sim_prefix) + ".fa",
            "--tree",
            str(sat_prefix) + ".sat.tree",
            "--model",
            MODEL,
            "--branch-id",
            branch_id,
            "--rate-file",
            str(sat_prefix) + ".rate",
            "--iqtree-report",
            str(sat_prefix) + ".iqtree",
            "--out",
            str(branch_reference_out),
            "--compare-sat-stat",
            sat_stat,
        ]
    )

    print(
        f"general_reference\tmodel=GTR+G4\ttaxa=6\tsplit=A,B\tbranch_id={branch_id}\tstatus=matched_native"
    )
    print(f"SatuTe generalized reference regression passed; outputs are in {outdir}")


if __name__ == "__main__":
    main()
