#!/usr/bin/env python3
"""Rebuild held-out null calibration from a merged truth audit stream."""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import tempfile
from collections import defaultdict
from contextlib import contextmanager
from pathlib import Path

from fdr_calibration_postprocess import (
    HELDOUT_FIELDS,
    calibration_cell,
    heldout_calibration,
    write_table,
)
from satute_analysis.fdr_calibration import TRUTH_BRANCH_FIELDS


@contextmanager
def open_truth_stream(path: Path):
    if path.suffix == ".zst":
        process = subprocess.Popen(
            ["zstd", "-dc", "--", str(path)],
            stdout=subprocess.PIPE,
            text=True,
            encoding="utf-8",
        )
        try:
            yield process.stdout
        finally:
            process.stdout.close()
            if process.wait() != 0:
                raise RuntimeError(f"zstd failed while reading {path}")
    else:
        with path.open("r", encoding="utf-8", newline="") as handle:
            yield handle


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth-audit", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--expected-training", type=int, default=2500)
    parser.add_argument("--expected-holdout", type=int, default=2500)
    args = parser.parse_args()
    if not 0.0 < args.alpha < 1.0:
        parser.error("--alpha must lie strictly between zero and one")

    training = defaultdict(list)
    holdout = defaultdict(list)
    with open_truth_stream(args.truth_audit) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(TRUTH_BRANCH_FIELDS):
            raise ValueError(f"Unexpected truth-audit schema: {reader.fieldnames}")
        for row in reader:
            if row["design"] != "exact_null" or row["truth"] != "null":
                continue
            target = training if row["sample"] == "training" else holdout
            target[calibration_cell(row)].append(float(row["p_value"]))

    rows = heldout_calibration(training, holdout, args.alpha)
    for row in rows:
        if row["training_nulls"] != args.expected_training:
            raise ValueError(f"Unexpected training-null count: {row}")
        if row["holdout_nulls"] != args.expected_holdout:
            raise ValueError(f"Unexpected holdout-null count: {row}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{args.out.name}.", suffix=".partial", dir=args.out.parent
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        write_table(temporary, HELDOUT_FIELDS, rows)
        os.replace(temporary, args.out)
    finally:
        temporary.unlink(missing_ok=True)
    print(f"Wrote {len(rows)} calibrated cells to {args.out}")


if __name__ == "__main__":
    main()
