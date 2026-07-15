"""Adapter and validators for native IQ-TREE SatuTe tabular output."""

from __future__ import annotations

import csv
import math
from pathlib import Path

from .benchmark_contracts import FORMULAS, SCHEMA_VERSION


REQUIRED_NATIVE_COLUMNS = {
    "ID",
    "Formula",
    "RateCategory",
    "LeftTaxa",
    "RightTaxa",
    "ValidSites",
    "SkippedSites",
    "satC",
    "satVar",
    "satSE",
    "satZ",
    "satP",
    "Alpha",
    "AlphaTaxonBonf",
    "Decision",
    "DecisionTaxonBonf",
    "FDR_BY",
    "DecisionFDR",
    "Length",
    "Split",
    "Modes",
    "Eigenvalues",
    "Weights",
}


def _probability(value: str) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) and 0.0 <= number <= 1.0 else None


def normalized_split(value) -> tuple[str, ...]:
    if isinstance(value, str):
        values = value.split(",")
    else:
        values = value
    return tuple(sorted(item.strip() for item in values if item and item.strip()))


def read_native_rows(path: str | Path, pooled_only: bool = True) -> list[dict]:
    """Read and strictly validate a native ``.sat.stat`` table."""

    path = Path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(
            csv.DictReader(
                (line for line in handle if line.strip() and not line.startswith("#")),
                delimiter="\t",
            )
        )
    if not rows:
        raise ValueError(f"Native SatuTe table contains no data rows: {path}")
    missing = REQUIRED_NATIVE_COLUMNS - set(rows[0])
    if missing:
        raise ValueError(f"Native SatuTe table is missing columns {sorted(missing)}: {path}")

    selected = [row for row in rows if not pooled_only or row["RateCategory"] == "pooled"]
    formulas = {row["Formula"] for row in selected}
    unexpected = formulas - set(FORMULAS)
    if unexpected:
        raise ValueError(f"Unexpected native formula families {sorted(unexpected)} in {path}")
    if pooled_only and formulas != set(FORMULAS):
        raise ValueError(f"Incomplete pooled formula families {sorted(formulas)} in {path}")

    for row in selected:
        for name in ("Decision", "DecisionTaxonBonf", "DecisionFDR"):
            if row[name] not in {"informative", "saturated", "not_tested"}:
                raise ValueError(f"Invalid {name}={row[name]!r} in {path}")
        if row["RateCategory"] == "pooled" and row["DecisionFDR"] == "not_tested":
            raise ValueError(f"Pooled row was not included in the FDR family in {path}")
    return selected


def index_target_rows(rows: list[dict], target_taxa) -> dict[str, dict]:
    target = normalized_split(target_taxa)
    matches = [row for row in rows if normalized_split(row.get("Split", "")) == target]
    indexed = {}
    for row in matches:
        formula = row["Formula"]
        if formula in indexed:
            raise ValueError(f"Duplicate pooled target row for formula {formula}: {target}")
        indexed[formula] = row
    return indexed


def target_detail(native: dict, alpha_fallback: float) -> dict:
    sat_p = _probability(native.get("satP"))
    left_taxa = int(native["LeftTaxa"])
    right_taxa = int(native["RightTaxa"])
    comparisons = left_taxa * right_taxa
    p_taxon_bonf = min(1.0, sat_p * comparisons) if sat_p is not None else ""
    return {
        "schema_version": SCHEMA_VERSION,
        "branch_id": native.get("ID", ""),
        "left_taxa": left_taxa,
        "right_taxa": right_taxa,
        "valid_sites": native.get("ValidSites", ""),
        "skipped_sites": native.get("SkippedSites", ""),
        "branch_length_used": native.get("Length", ""),
        "alpha": native.get("Alpha", alpha_fallback),
        "satC": native.get("satC", ""),
        "satVar": native.get("satVar", ""),
        "satSE": native.get("satSE", ""),
        "satZ": native.get("satZ", ""),
        "satP": native.get("satP", ""),
        "decision_unadjusted": native.get("Decision", ""),
        "alpha_taxon_bonf": native.get("AlphaTaxonBonf", ""),
        "p_taxon_bonf": p_taxon_bonf,
        "decision_taxon_bonf": native.get("DecisionTaxonBonf", ""),
        "fdr_by": native.get("FDR_BY", ""),
        "decision_fdr": native.get("DecisionFDR", ""),
        "modes": native.get("Modes", ""),
        "eigenvalues": native.get("Eigenvalues", ""),
        "weights": native.get("Weights", ""),
    }


def fdr_family_audit(rows: list[dict], formula: str) -> dict:
    family = [row for row in rows if row["Formula"] == formula]
    if not family:
        raise ValueError(f"No pooled native rows for formula {formula}")
    branch_ids = [row["ID"] for row in family]
    if len(branch_ids) != len(set(branch_ids)):
        raise ValueError(f"Duplicate branch IDs in {formula} FDR family")
    p_values = [_probability(row.get("satP")) for row in family]
    adjusted = [_probability(row.get("FDR_BY")) for row in family]
    return {
        "schema_version": SCHEMA_VERSION,
        "formula": formula,
        "family_size": len(family),
        "valid_p_values": sum(value is not None for value in p_values),
        "unadjusted_rejections": sum(row["Decision"] == "informative" for row in family),
        "by_fdr_rejections": sum(row["DecisionFDR"] == "informative" for row in family),
        "minimum_p_value": min((value for value in p_values if value is not None), default=""),
        "minimum_fdr_by": min((value for value in adjusted if value is not None), default=""),
    }


def branch_audit_detail(rows: list[dict], target_taxa) -> list[dict]:
    """Return compact pooled branch rows suitable for empirical FDR audits."""

    target = normalized_split(target_taxa)
    detail = []
    for row in rows:
        sat_p = _probability(row.get("satP"))
        left_taxa = int(row["LeftTaxa"])
        right_taxa = int(row["RightTaxa"])
        comparisons = left_taxa * right_taxa
        detail.append(
            {
                "branch_id": row.get("ID", ""),
                "split": row.get("Split", ""),
                "left_taxa": left_taxa,
                "right_taxa": right_taxa,
                "satP": row.get("satP", ""),
                "decision_unadjusted": row.get("Decision", ""),
                "alpha_taxon_bonf": row.get("AlphaTaxonBonf", ""),
                "p_taxon_bonf": min(1.0, sat_p * comparisons) if sat_p is not None else "",
                "decision_taxon_bonf": row.get("DecisionTaxonBonf", ""),
                "fdr_by": row.get("FDR_BY", ""),
                "decision_fdr": row.get("DecisionFDR", ""),
                "is_target": int(normalized_split(row.get("Split", "")) == target),
            }
        )
    return detail
