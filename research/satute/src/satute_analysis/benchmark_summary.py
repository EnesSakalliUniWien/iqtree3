"""Streaming aggregation and validation for benchmark TSV artifacts."""

from __future__ import annotations

import csv
import os
from collections import defaultdict
from pathlib import Path

from .benchmark_contracts import (
    BASE_DETAIL_FIELDS,
    SCHEMA_VERSION,
    SUMMARY_FIELDS,
    applicable_decision_rules,
)


def validate_tsv_header(path: str | Path, expected_fields=BASE_DETAIL_FIELDS) -> None:
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return
    with path.open("r", encoding="utf-8", newline="") as handle:
        observed = next(csv.reader(handle, delimiter="\t"), [])
    if tuple(observed) != tuple(expected_fields):
        raise ValueError(
            f"Refusing to resume incompatible TSV schema in {path}; "
            f"expected version {SCHEMA_VERSION} fields"
        )


def aggregate_detail(detail_path: str | Path, summary_path: str | Path) -> None:
    """Create a normalized decision-rule summary using an atomic replacement."""

    groups = defaultdict(lambda: {"evaluated": 0, "informative": 0, "missing": 0})
    with Path(detail_path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(BASE_DETAIL_FIELDS):
            raise ValueError(f"Unexpected benchmark detail schema: {detail_path}")
        for row in reader:
            if int(row["schema_version"]) != SCHEMA_VERSION:
                raise ValueError(f"Mixed schema versions in {detail_path}")
            for rule in applicable_decision_rules(row["scenario"]):
                key = (
                    row["tree_case"],
                    row["simulation_model"],
                    row["evaluation_model"],
                    row["nsites"],
                    row["branch_length"],
                    row["scenario"],
                    rule.name,
                    row["formula"],
                )
                if row["target_found"] == "1":
                    groups[key]["evaluated"] += 1
                    groups[key]["informative"] += row[rule.decision_column] == "informative"
                else:
                    groups[key]["missing"] += 1

    summary_path = Path(summary_path)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = summary_path.with_name(f".{summary_path.name}.tmp.{os.getpid()}")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS, delimiter="\t")
            writer.writeheader()
            for key in sorted(
                groups,
                key=lambda value: (
                    value[0], value[1], value[2], int(value[3]), float(value[4]),
                    value[5], value[6], value[7],
                ),
            ):
                group = groups[key]
                evaluated = group["evaluated"]
                writer.writerow(
                    {
                        "schema_version": SCHEMA_VERSION,
                        "tree_case": key[0],
                        "simulation_model": key[1],
                        "evaluation_model": key[2],
                        "nsites": key[3],
                        "branch_length": key[4],
                        "scenario": key[5],
                        "decision_rule": key[6],
                        "formula": key[7],
                        "evaluated": evaluated,
                        "informative": group["informative"],
                        "missing_split": group["missing"],
                        "fraction_informative": group["informative"] / evaluated if evaluated else "",
                    }
                )
        os.replace(temporary, summary_path)
    finally:
        if temporary.exists():
            temporary.unlink()
