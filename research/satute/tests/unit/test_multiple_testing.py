#!/usr/bin/env python3

import math

from satute_analysis.multiple_testing import (
    annotate_satute_rows,
    benjamini_yekutieli,
    paper_taxon_bonferroni,
)


def assert_close(observed, expected, tolerance=1e-14):
    if not math.isclose(observed, expected, rel_tol=0.0, abs_tol=tolerance):
        raise AssertionError(f"observed {observed} != expected {expected}")


def check_paper_bonferroni():
    result = paper_taxon_bonferroni(0.001, 8, 8, 0.05)
    assert_close(result["AlphaTaxonBonf"], 0.05 / 64.0)
    assert_close(result["PTaxonBonf"], 0.064)
    if result["DecisionTaxonBonf"] != "saturated":
        raise AssertionError(result)


def check_by_adjustment():
    adjusted = benjamini_yekutieli([0.001, 0.01, 0.2])
    expected = [0.0055, 0.0275, 0.3666666666666667]
    for observed, target in zip(adjusted, expected):
        assert_close(observed, target)

    with_missing = benjamini_yekutieli([0.001, None, 0.01])
    assert_close(with_missing[0], 0.0055)
    if with_missing[1] is not None:
        raise AssertionError(with_missing)
    assert_close(with_missing[2], 0.0275)


def check_separate_formula_families():
    rows = [
        {
            "Split": split,
            "Formula": formula,
            "RateCategory": category,
            "satP": p_value,
            "LeftTaxa": left_taxa,
            "RightTaxa": right_taxa,
        }
        for split, formula, category, p_value, left_taxa, right_taxa in (
            ("A", "dominant", "pooled", 0.001, 1, 3),
            ("C,D", "dominant", "pooled", 0.04, 2, 2),
            ("A", "eigenvalue_weighted", "pooled", 0.02, 1, 3),
            ("C,D", "eigenvalue_weighted", "pooled", 0.03, 2, 2),
            ("A", "dominant", "1", 0.0001, 1, 3),
            ("A", "eigenvalue_weighted", "1", 0.0002, 1, 3),
        )
    ]
    original_rows = [dict(row) for row in rows]
    annotated = annotate_satute_rows(rows, 0.05)
    if rows != original_rows:
        raise AssertionError("annotation mutated its input rows")

    dominant = [row for row in annotated if row["Formula"] == "dominant" and row["RateCategory"] == "pooled"]
    weighted = [
        row
        for row in annotated
        if row["Formula"] == "eigenvalue_weighted" and row["RateCategory"] == "pooled"
    ]
    dominant_expected = benjamini_yekutieli([0.001, 0.04])
    weighted_expected = benjamini_yekutieli([0.02, 0.03])
    for row, expected in zip(dominant, dominant_expected):
        assert_close(row["FDR_BY"], expected)
    for row, expected in zip(weighted, weighted_expected):
        assert_close(row["FDR_BY"], expected)

    categories = [row for row in annotated if row["RateCategory"] != "pooled"]
    if any(row["FDR_BY"] is not None or row["DecisionFDR"] != "not_tested" for row in categories):
        raise AssertionError(categories)
    if any(row["PTaxonBonf"] is None for row in annotated):
        raise AssertionError("paper correction was not applied to every valid row")


def main():
    check_paper_bonferroni()
    check_by_adjustment()
    check_separate_formula_families()
    print("Multiple-testing unit checks passed")


if __name__ == "__main__":
    main()
