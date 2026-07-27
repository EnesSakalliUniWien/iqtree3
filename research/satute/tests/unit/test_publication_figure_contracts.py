#!/usr/bin/env python3
"""Regression checks for the scientific contracts of publication figures."""

from satute_analysis.benchmark_contracts import applicable_decision_rules
from satute_analysis.publication_figures import (
    CROSS_MODEL_CONTRACT,
    MISSPECIFICATION_CONTRACT,
    stable_case_summary,
    stable_figure_root,
)


def allowed_rules(scenario: str) -> set[str]:
    return {rule.name for rule in applicable_decision_rules(scenario)}


def main() -> None:
    for series in MISSPECIFICATION_CONTRACT.series:
        if series.decision_rule not in allowed_rules(series.scenario):
            raise SystemExit(
                f"Figure 5 uses invalid scenario/rule pair: "
                f"{series.scenario}/{series.decision_rule}"
            )

    if CROSS_MODEL_CONTRACT.decision_rule not in allowed_rules(
        CROSS_MODEL_CONTRACT.scenario
    ):
        raise SystemExit(
            "Figure 6 uses an invalid scenario/decision-rule combination"
        )
    if CROSS_MODEL_CONTRACT.site_lengths != (100, 1000):
        raise SystemExit("Figure 6 must compare the common 100- and 1,000-site grids")

    for case_id in (
        "jc_JC",
        "gtr_GTR_PF06346_G4",
        "protein_LG_G4",
        "misspec_JC",
    ):
        path = stable_case_summary(case_id)
        if "artifacts/publication/power" not in path.as_posix():
            raise SystemExit(f"Figure path is not case-local: {path}")
        if "artifacts/cluster" in path.as_posix():
            raise SystemExit(f"Figure path exposes a hashed cluster release: {path}")

    for contract in (MISSPECIFICATION_CONTRACT, CROSS_MODEL_CONTRACT):
        root = stable_figure_root(contract.id)
        if root.parts[-4:-1] != ("artifacts", "publication", "figures"):
            raise SystemExit(f"Figure output is not publication-local: {root}")

    print("Publication figure scientific contracts: passed")


if __name__ == "__main__":
    main()
