#!/usr/bin/env python3
"""Summarize per-replicate timing and identify benchmark bottlenecks."""

from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path


def percentile(values: list[float], probability: float) -> float:
    if not values:
        return float("nan")
    ordered = sorted(values)
    index = min(len(ordered) - 1, int(round(probability * (len(ordered) - 1))))
    return ordered[index]


def main() -> None:
    parser = argparse.ArgumentParser(description="Profile a merged SatuTe timing table.")
    parser.add_argument("--timing", required=True)
    args = parser.parse_args()

    rows = list(csv.DictReader(Path(args.timing).open(encoding="utf-8"), delimiter="\t"))
    if not rows:
        raise SystemExit("Timing table is empty")
    metrics = ["simulation_seconds", "true_fixed_seconds", "true_ml_lengths_seconds", "ml_tree_seconds", "total_seconds"]
    print(f"replicates={len(rows)}")
    cache_hits = sum(row.get("simulation_cache_hit") == "1" for row in rows)
    print(f"simulation_cache_hit_fraction={cache_hits / len(rows):.4f}")
    for metric in metrics:
        values = [float(row[metric]) for row in rows if row.get(metric, "") not in {"", "NA"}]
        if not values:
            continue
        print(
            f"{metric}: median={statistics.median(values):.4f}s "
            f"p95={percentile(values, 0.95):.4f}s "
            f"mean={statistics.mean(values):.4f}s"
        )
    total = sum(float(row["total_seconds"]) for row in rows)
    print(f"stage_share_simulation={sum(float(row['simulation_seconds']) for row in rows) / total:.4f}")
    print(f"stage_share_true_fixed={sum(float(row['true_fixed_seconds']) for row in rows) / total:.4f}")
    print(f"stage_share_true_ml_lengths={sum(float(row['true_ml_lengths_seconds']) for row in rows) / total:.4f}")
    print(f"stage_share_ml_tree={sum(float(row['ml_tree_seconds']) for row in rows) / total:.4f}")


if __name__ == "__main__":
    main()
