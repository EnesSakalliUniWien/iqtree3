"""Multiple-testing corrections shared by SatuTe Python workflows.

The paper's taxon-pair Bonferroni correction and the tree-wide
Benjamini--Yekutieli correction address different hypothesis families.  They
are therefore calculated independently and reported side by side.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping


SATUTE_FORMULAS = ("dominant", "eigenvalue_weighted")


def _valid_probability(value) -> bool:
    return isinstance(value, (int, float)) and math.isfinite(value) and 0.0 <= value <= 1.0


def paper_taxon_bonferroni(p_value, left_taxa: int, right_taxa: int, alpha: float) -> dict:
    """Return the paper's ``alpha / (n_left * n_right)`` correction.

    The returned adjusted p-value is the equivalent Bonferroni representation
    of the corrected-alpha decision.  Undefined p-values retain the planned
    threshold but receive no adjusted p-value or decision.
    """

    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must be between 0 and 1")
    if left_taxa < 1 or right_taxa < 1:
        raise ValueError("both branch sides must contain at least one taxon")

    comparisons = left_taxa * right_taxa
    adjusted_alpha = alpha / comparisons
    if not _valid_probability(p_value):
        return {
            "AlphaTaxonBonf": adjusted_alpha,
            "PTaxonBonf": None,
            "DecisionTaxonBonf": "not_tested",
        }

    adjusted_p = min(1.0, float(p_value) * comparisons)
    return {
        "AlphaTaxonBonf": adjusted_alpha,
        "PTaxonBonf": adjusted_p,
        "DecisionTaxonBonf": "informative" if p_value <= adjusted_alpha else "saturated",
    }


def benjamini_yekutieli(p_values: Iterable[float | None]) -> list[float | None]:
    """Return monotone BY-adjusted p-values in the original input order.

    Every supplied position belongs to the planned family.  Undefined p-values
    are ranked as one so that they cannot be rejected or shrink the
    multiplicity penalty, but their returned adjusted value remains ``None``.
    Stable input indices resolve ties deterministically.
    """

    values = list(p_values)
    family_size = len(values)
    if family_size == 0:
        return []

    ordered = sorted(
        (
            (float(value) if _valid_probability(value) else 1.0, index)
            for index, value in enumerate(values)
        ),
        key=lambda item: (item[0], item[1]),
    )
    harmonic = sum(1.0 / rank for rank in range(1, family_size + 1))
    adjusted = [None] * family_size
    running_minimum = 1.0

    for offset in range(family_size - 1, -1, -1):
        rank = offset + 1
        p_value, original_index = ordered[offset]
        candidate = min(1.0, p_value * family_size * harmonic / rank)
        running_minimum = min(running_minimum, candidate)
        if _valid_probability(values[original_index]):
            adjusted[original_index] = running_minimum

    return adjusted


def annotate_satute_rows(
    rows: Iterable[Mapping],
    alpha: float,
    formulas: Iterable[str] = SATUTE_FORMULAS,
) -> list[dict]:
    """Add paper Bonferroni and separate formula-specific BY results.

    Required row fields are ``Formula``, ``RateCategory``, ``satP``,
    ``LeftTaxa`` and ``RightTaxa``.  The input rows are copied.  The BY family
    for each formula contains only its pooled branch rows; rate-category rows
    receive ``DecisionFDR=not_tested``.
    """

    formula_order = tuple(formulas)
    if len(formula_order) != len(set(formula_order)):
        raise ValueError("formula families must be unique")

    annotated = []
    for source_row in rows:
        row = dict(source_row)
        formula = row.get("Formula")
        if formula not in formula_order:
            raise ValueError(f"unsupported SatuTe formula family: {formula}")
        row["Alpha"] = alpha
        row.update(
            paper_taxon_bonferroni(
                row.get("satP"),
                int(row["LeftTaxa"]),
                int(row["RightTaxa"]),
                alpha,
            )
        )
        row["FDR_BY"] = None
        row["DecisionFDR"] = "not_tested"
        annotated.append(row)

    for formula in formula_order:
        family_indices = [
            index
            for index, row in enumerate(annotated)
            if row["Formula"] == formula and row["RateCategory"] == "pooled"
        ]
        adjusted = benjamini_yekutieli(annotated[index].get("satP") for index in family_indices)
        for row_index, adjusted_p in zip(family_indices, adjusted):
            row = annotated[row_index]
            row["FDR_BY"] = adjusted_p
            if adjusted_p is not None:
                row["DecisionFDR"] = "informative" if adjusted_p <= alpha else "saturated"

    return annotated
