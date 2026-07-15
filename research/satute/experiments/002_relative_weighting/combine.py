#!/usr/bin/env python3
"""Combine and validate schema-versioned SatuTe shard artifacts."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

from satute_analysis.benchmark_contracts import (
    BASE_DETAIL_FIELDS,
    BRANCH_AUDIT_FIELDS,
    FDR_AUDIT_FIELDS,
    TIMING_FIELDS,
)
from satute_analysis.benchmark_summary import aggregate_detail


TABLES = {
    "head_to_head_detail.tsv": BASE_DETAIL_FIELDS,
    "head_to_head_branch_audit.tsv": BRANCH_AUDIT_FIELDS,
    "head_to_head_fdr_audit.tsv": FDR_AUDIT_FIELDS,
    "head_to_head_timing.tsv": TIMING_FIELDS,
}

KEY_FIELDS = {
    "head_to_head_detail.tsv": (
        "schema_version", "tree_case", "simulation_model", "evaluation_model", "nsites",
        "branch_length", "replicate", "scenario", "formula",
    ),
    "head_to_head_branch_audit.tsv": (
        "schema_version", "tree_case", "simulation_model", "evaluation_model", "nsites",
        "branch_length", "replicate", "scenario", "formula", "branch_id",
    ),
    "head_to_head_fdr_audit.tsv": (
        "schema_version", "tree_case", "simulation_model", "evaluation_model", "nsites",
        "branch_length", "replicate", "scenario", "formula",
    ),
    "head_to_head_timing.tsv": (
        "schema_version", "tree_case", "simulation_model", "evaluation_model", "nsites",
        "branch_length", "replicate",
    ),
}


def combine_table(inputs: list[Path], output: Path, expected_fields: tuple[str, ...], key_fields: tuple[str, ...]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp.{os.getpid()}")
    seen: set[tuple[str, ...]] = set()
    try:
        with temporary.open("w", encoding="utf-8", newline="") as out_handle:
            writer = None
            for input_path in inputs:
                with input_path.open("r", encoding="utf-8", newline="") as in_handle:
                    reader = csv.DictReader(in_handle, delimiter="\t")
                    observed = tuple(reader.fieldnames or ())
                    if observed != tuple(expected_fields):
                        raise ValueError(
                            f"Incompatible schema in {input_path}: expected {expected_fields}, got {observed}"
                        )
                    if writer is None:
                        writer = csv.DictWriter(out_handle, fieldnames=expected_fields, delimiter="\t")
                        writer.writeheader()
                    for row in reader:
                        if "replicate" in row:
                            logical_key = tuple(row.get(name, "") for name in key_fields)
                            if logical_key in seen:
                                raise ValueError(f"Duplicate logical row while combining: {input_path}: {logical_key}")
                            seen.add(logical_key)
                        writer.writerow(row)
        os.replace(temporary, output)
    finally:
        if temporary.exists():
            temporary.unlink()


def main() -> None:
    parser = argparse.ArgumentParser(description="Combine schema-versioned SatuTe shard artifacts.")
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    detail_inputs = [Path(value) for value in args.inputs]
    for path in detail_inputs:
        if not path.is_file():
            raise SystemExit(f"Missing detail shard: {path}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for name, fields in TABLES.items():
        auxiliary = [path.parent / name for path in detail_inputs]
        missing = [path for path in auxiliary if not path.is_file()]
        if missing:
            raise SystemExit(f"Missing companion {name} files: {missing[:3]}")
        combine_table(auxiliary, outdir / name, fields, KEY_FIELDS[name])

    aggregate_detail(outdir / "head_to_head_detail.tsv", outdir / "head_to_head_summary.tsv")
    print(f"Detail:       {outdir / 'head_to_head_detail.tsv'}")
    print(f"Summary:      {outdir / 'head_to_head_summary.tsv'}")
    print(f"Branch audit: {outdir / 'head_to_head_branch_audit.tsv'}")
    print(f"FDR audit:    {outdir / 'head_to_head_fdr_audit.tsv'}")
    print(f"Timing:       {outdir / 'head_to_head_timing.tsv'}")


if __name__ == "__main__":
    main()
