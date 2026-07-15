#!/usr/bin/env python3

import csv
import tempfile
from pathlib import Path

from satute_analysis.benchmark_contracts import (
    BASE_DETAIL_FIELDS,
    BRANCH_AUDIT_FIELDS,
    FDR_AUDIT_FIELDS,
    SUMMARY_FIELDS,
    TIMING_FIELDS,
    applicable_decision_rules,
)
from satute_analysis.benchmark_summary import aggregate_detail


def make_detail(path: Path) -> None:
    row = {field: "" for field in BASE_DETAIL_FIELDS}
    row.update(
        {
            "schema_version": "2",
            "tree_case": "five_external",
            "simulation_model": "JC",
            "evaluation_model": "JC",
            "fitted_evaluation_model": "JC",
            "nsites": "100",
            "branch_length": "0.5",
            "replicate": "1",
            "seed": "7",
            "simulator": "alisim",
            "target_split": "B",
            "scenario": "true_tree_fixed_lengths",
            "implementation": "iqtree_native",
            "target_found": "1",
            "satP": "0.01",
            "decision_unadjusted": "informative",
            "alpha_taxon_bonf": "0.0125",
            "p_taxon_bonf": "0.04",
            "decision_taxon_bonf": "saturated",
            "fdr_by": "0.03",
            "decision_fdr": "informative",
        }
    )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=BASE_DETAIL_FIELDS, delimiter="\t")
        writer.writeheader()
        for formula in ("dominant", "eigenvalue_weighted"):
            row["formula"] = formula
            writer.writerow(row)


def main() -> None:
    schemas = [BASE_DETAIL_FIELDS, BRANCH_AUDIT_FIELDS, FDR_AUDIT_FIELDS, SUMMARY_FIELDS, TIMING_FIELDS]
    if any(len(fields) != len(set(fields)) for fields in schemas):
        raise SystemExit("Schema field names must be unique within each table")
    if [rule.name for rule in applicable_decision_rules("ml_tree")] != [
        "unadjusted", "taxon_bonferroni", "by_fdr"
    ]:
        raise SystemExit("ML-tree decision-rule contract changed unexpectedly")
    if [rule.name for rule in applicable_decision_rules("true_tree_fixed_lengths")] != [
        "unadjusted", "by_fdr"
    ]:
        raise SystemExit("Fixed-tree decision-rule contract changed unexpectedly")

    with tempfile.TemporaryDirectory() as temporary:
        detail = Path(temporary) / "detail.tsv"
        summary = Path(temporary) / "summary.tsv"
        make_detail(detail)
        aggregate_detail(detail, summary)
        with summary.open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        observed = {(row["scenario"], row["decision_rule"], row["formula"]) for row in rows}
        expected = {
            ("true_tree_fixed_lengths", rule, formula)
            for rule in ("unadjusted", "by_fdr")
            for formula in ("dominant", "eigenvalue_weighted")
        }
        if observed != expected:
            raise SystemExit(f"Unexpected normalized summary keys: {sorted(observed)}")
    print("Benchmark schema contracts: passed")


if __name__ == "__main__":
    main()
