#!/usr/bin/env python3
"""Verify the SatuTe saturated-null expectation by exact pattern enumeration.

Requires NumPy.  The script enumerates all 4^4 patterns for the fixed
four-taxon GTR+G4 diagnostic model, proves the hard and soft score means are
zero numerically, then checks their Studentized null distributions by direct
multinomial sampling from that exact distribution.
"""

import argparse
import csv
import itertools
import math
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from satute_analysis.models import build_q, reversible_eigendecomposition, transition_matrix


MODEL = "GTR{1,2,1,1,2,1}+F{0.30,0.20,0.20,0.30}"
RATES = np.array([0.03338775338, 0.2519159176, 0.8202684819, 2.894427847])
PROPORTIONS = np.repeat(0.25, 4)


def exact_patterns():
    q, pi = build_q(MODEL)
    eigenvalues, eigenvectors, inverse = reversible_eigendecomposition(q, pi)
    transitions = [
        transition_matrix(eigenvalues, eigenvectors, inverse, 0.05 * rate)
        for rate in RATES
    ]
    zero = int(np.argmin(np.abs(eigenvalues)))
    dominant = max((index for index in range(4) if index != zero), key=lambda index: eigenvalues[index])

    probabilities = []
    hard_category = []
    hard_left = []
    hard_right = []
    soft_score = []
    for pattern in itertools.product(range(4), repeat=4):
        bases, left_factors, right_factors = [], [], []
        for category in range(4):
            left = pi * transitions[category][:, pattern[0]] * transitions[category][:, pattern[1]]
            right = pi * transitions[category][:, pattern[2]] * transitions[category][:, pattern[3]]
            left_base, right_base = left.sum(), right.sum()
            bases.append(PROPORTIONS[category] * left_base * right_base)
            left_factors.append(float(eigenvectors[:, dominant] @ (left / left_base)))
            right_factors.append(float(eigenvectors[:, dominant] @ (right / right_base)))
        probability = sum(bases)
        responsibilities = np.array(bases) / probability
        assigned = int(np.argmax(responsibilities))
        probabilities.append(probability)
        hard_category.append(assigned)
        hard_left.append(left_factors[assigned])
        hard_right.append(right_factors[assigned])
        # At the saturated limit the slowest category is the only absolute
        # component that remains after division by its common tiny weight.
        soft_score.append(responsibilities[0] * left_factors[0] * right_factors[0])

    probabilities = np.array(probabilities)
    probabilities /= probabilities.sum()
    hard_category = np.array(hard_category)
    hard_left = np.array(hard_left)
    hard_right = np.array(hard_right)
    return probabilities, hard_category, hard_left, hard_right, np.array(soft_score)


def simulate(probabilities, categories, left, right, soft, nsites, replicates, rng):
    hard_score = left * right
    hard_z, soft_z = [], []
    for start in range(0, replicates, 1000):
        count = min(1000, replicates - start)
        samples = rng.multinomial(nsites, probabilities, size=count)
        hard_mean = samples @ hard_score / nsites
        hard_variance = np.zeros(count)
        for category in sorted(set(categories)):
            mask = categories == category
            category_n = samples[:, mask].sum(axis=1)
            valid = category_n > 0
            left_second = np.zeros(count)
            right_second = np.zeros(count)
            left_second[valid] = samples[valid][:, mask] @ (left[mask] ** 2) / category_n[valid]
            right_second[valid] = samples[valid][:, mask] @ (right[mask] ** 2) / category_n[valid]
            hard_variance += category_n / nsites * left_second * right_second
        hard_z.extend(hard_mean / np.sqrt(hard_variance / nsites))

        soft_mean = samples @ soft / nsites
        soft_variance = (samples @ (soft * soft) - nsites * soft_mean * soft_mean) / (nsites - 1)
        soft_z.extend(soft_mean / np.sqrt(soft_variance / nsites))
    return np.array(hard_z), np.array(soft_z)


def summary_row(nsites, method, values):
    return {
        "nsites": nsites,
        "method": method,
        "replicates": len(values),
        "mean_z": float(values.mean()),
        "sd_z": float(values.std(ddof=1)),
        "q95_z": float(np.quantile(values, 0.95)),
        "nominal_rejection": float(np.mean(values > 1.6448536269514722)),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--site-lengths", default="100,500,1000,5000,10000")
    parser.add_argument("--replicates", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=20260710)
    parser.add_argument(
        "--output",
        default=str(PROJECT_ROOT / "artifacts/local/exact-null-verification.tsv"),
    )
    args = parser.parse_args()

    probabilities, categories, left, right, soft = exact_patterns()
    hard = left * right
    exact_hard_mean = float(probabilities @ hard)
    exact_soft_mean = float(probabilities @ soft)
    if abs(probabilities.sum() - 1.0) > 1e-12:
        raise AssertionError("Pattern probabilities do not sum to one")
    if abs(exact_hard_mean) > 1e-12 or abs(exact_soft_mean) > 1e-12:
        raise AssertionError(f"Nonzero exact null mean: hard={exact_hard_mean}, soft={exact_soft_mean}")

    rng = np.random.default_rng(args.seed)
    rows = []
    for nsites in [int(value) for value in args.site_lengths.split(",") if value.strip()]:
        hard_z, soft_z = simulate(probabilities, categories, left, right, soft, nsites, args.replicates, rng)
        rows.extend((summary_row(nsites, "hard", hard_z), summary_row(nsites, "soft", soft_z)))

    for row in rows:
        if abs(row["mean_z"]) > 0.05:
            raise AssertionError(f"Null centering failed: {row}")
        if not 0.95 <= row["sd_z"] <= 1.05:
            raise AssertionError(f"Null scale failed: {row}")
        if not 0.04 <= row["nominal_rejection"] <= 0.06:
            raise AssertionError(f"Null rejection failed: {row}")

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"exact hard mean: {exact_hard_mean:.3g}")
    print(f"exact soft mean: {exact_soft_mean:.3g}")
    print(output)


if __name__ == "__main__":
    main()
