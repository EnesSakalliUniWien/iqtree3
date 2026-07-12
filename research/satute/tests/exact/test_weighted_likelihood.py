#!/usr/bin/env python3
"""Exact checks for the refined weighted SatuTe likelihood construction."""

import itertools
import math
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from satute_analysis.models import build_q, reversible_eigendecomposition, transition_matrix


MODEL = "GTR{1,2,1,1,2,1}+F{0.30,0.20,0.20,0.30}"
VARIABLE_RATES = np.array([0.03338775338, 0.2519159176, 0.8202684819, 2.894427847])
PATTERNS = tuple(itertools.product(range(4), repeat=4))


def category_terms(eigenvalues, eigenvectors, inverse, pi, modes, rate):
    side_transition = transition_matrix(eigenvalues, eigenvectors, inverse, 0.05 * rate)
    left_partials = []
    right_partials = []
    independent = []
    coherence = []
    invariant_joint = []

    for pattern in PATTERNS:
        left = side_transition[:, pattern[0]] * side_transition[:, pattern[1]]
        right = side_transition[:, pattern[2]] * side_transition[:, pattern[3]]
        left_probability = float(pi @ left)
        right_probability = float(pi @ right)
        left_partials.append(left)
        right_partials.append(right)
        independent.append(left_probability * right_probability)
        invariant_joint.append(float(pi @ (left * right)))
        if left_probability > 0.0 and right_probability > 0.0:
            left_posterior = pi * left / left_probability
            right_posterior = pi * right / right_probability
            left_factors = eigenvectors[:, modes].T @ left_posterior
            right_factors = eigenvectors[:, modes].T @ right_posterior
            coherence.append(left_factors * right_factors)
        else:
            coherence.append(np.zeros(len(modes)))

    return {
        "left": np.asarray(left_partials),
        "right": np.asarray(right_partials),
        "independent": np.asarray(independent),
        "coherence": np.asarray(coherence),
        "invariant_joint": np.asarray(invariant_joint),
    }


def check_mixture(eigenvalues, eigenvectors, inverse, pi, modes, rates, proportions, branch_length):
    categories = [
        category_terms(eigenvalues, eigenvectors, inverse, pi, modes, float(rate))
        for rate in rates
    ]
    positive_exponents = [
        eigenvalues[mode] * branch_length * rate
        for rate, proportion in zip(rates, proportions)
        if proportion > 0.0 and rate > 0.0
        for mode in modes
    ]
    shift = max(positive_exponents, default=0.0)
    null_probability = np.zeros(len(PATTERNS))
    finite_probability = np.zeros(len(PATTERNS))
    scaled_numerator = np.zeros(len(PATTERNS))

    for rate, proportion, category in zip(rates, proportions, categories):
        if rate == 0.0:
            null_probability += proportion * category["invariant_joint"]
        else:
            null_probability += proportion * category["independent"]
            scaled_weights = np.exp(eigenvalues[modes] * branch_length * rate - shift)
            scaled_numerator += (
                proportion
                * category["independent"]
                * (category["coherence"] @ scaled_weights)
            )

        focal_transition = transition_matrix(
            eigenvalues,
            eigenvectors,
            inverse,
            branch_length * rate,
        )
        finite_probability += proportion * np.asarray(
            [
                np.einsum("a,a,ab,b->", pi, left, focal_transition, right)
                for left, right in zip(category["left"], category["right"])
            ]
        )

    if np.any(null_probability <= 0.0):
        raise AssertionError("The exact mixture null contains a zero-probability pattern")
    null_probability /= null_probability.sum()
    finite_probability /= finite_probability.sum()
    scaled_score = scaled_numerator / null_probability
    ratio_error = float(
        np.max(
            np.abs(
                finite_probability / null_probability
                - 1.0
                - math.exp(shift) * scaled_score
            )
        )
    )
    null_mean = float(null_probability @ scaled_score)
    if ratio_error > 2e-11:
        raise AssertionError(f"Weighted likelihood-ratio identity failed: {ratio_error}")
    if abs(null_mean) > 2e-11:
        raise AssertionError(f"Weighted score null centering failed: {null_mean}")
    return shift, ratio_error, null_mean


def check_linear_optimum(eigenvalues, eigenvectors, inverse, pi, modes):
    category = category_terms(eigenvalues, eigenvectors, inverse, pi, modes, 1.0)
    null_probability = category["independent"]
    null_probability /= null_probability.sum()
    coherence = category["coherence"]
    covariance = np.einsum("s,si,sj->ij", null_probability, coherence, coherence)
    branch_length = 0.5
    persistence = np.exp(eigenvalues[modes] * branch_length)
    focal_transition = transition_matrix(
        eigenvalues, eigenvectors, inverse, branch_length
    )
    finite_probability = np.asarray(
        [
            np.einsum("a,a,ab,b->", pi, left, focal_transition, right)
            for left, right in zip(category["left"], category["right"])
        ]
    )
    finite_probability /= finite_probability.sum()
    alternative_mean = finite_probability @ coherence
    mean_error = float(np.max(np.abs(alternative_mean - covariance @ persistence)))
    if mean_error > 2e-12:
        raise AssertionError(f"Alternative mean identity failed: {mean_error}")

    gls = np.linalg.solve(covariance + 1e-12 * np.eye(len(modes)), persistence)

    def noncentrality(weights):
        return float(
            weights @ alternative_mean
            / math.sqrt(weights @ covariance @ weights)
        )

    persistence_ncp = noncentrality(persistence)
    gls_ncp = noncentrality(gls)
    if persistence_ncp + 1e-10 < gls_ncp:
        raise AssertionError("Covariance-whitened weights exceeded the exact persistence optimum")
    return mean_error, persistence_ncp, gls_ncp


def main():
    q, pi = build_q(MODEL)
    eigenvalues, eigenvectors, inverse = reversible_eigendecomposition(q, pi)
    zero = int(np.argmin(np.abs(eigenvalues)))
    modes = [index for index in range(len(eigenvalues)) if index != zero]

    mean_error, persistence_ncp, gls_ncp = check_linear_optimum(
        eigenvalues, eigenvectors, inverse, pi, modes
    )
    print(
        "linear optimum: "
        f"mean error={mean_error:.3g}; persistence NCP={persistence_ncp:.8g}; "
        f"GLS NCP={gls_ncp:.8g}"
    )

    for branch_length in (0.5, 4.0, 50.0):
        shift, ratio_error, null_mean = check_mixture(
            eigenvalues,
            eigenvectors,
            inverse,
            pi,
            modes,
            VARIABLE_RATES,
            np.repeat(0.25, 4),
            branch_length,
        )
        print(
            f"positive-rate mixture t={branch_length:g}: shift={shift:.8g}; "
            f"ratio error={ratio_error:.3g}; null mean={null_mean:.3g}"
        )

    invariant_rates = np.concatenate(([0.0], VARIABLE_RATES))
    invariant_proportions = np.concatenate(([0.1], np.repeat(0.225, 4)))
    for branch_length in (0.5, 4.0, 50.0, 10000.0):
        shift, ratio_error, null_mean = check_mixture(
            eigenvalues,
            eigenvectors,
            inverse,
            pi,
            modes,
            invariant_rates,
            invariant_proportions,
            branch_length,
        )
        print(
            f"+I mixture t={branch_length:g}: shift={shift:.8g}; "
            f"ratio error={ratio_error:.3g}; null mean={null_mean:.3g}"
        )


if __name__ == "__main__":
    main()
