#!/usr/bin/env python3
"""Exact checks for the refined weighted SatuTe likelihood construction.

Besides the likelihood-ratio identities, this verifier checks the local
geometry of the independence null.  Modal coherence scores must lie in the
orthogonal complement of the nuisance tangent space, their Gram matrix must
define non-degenerate local information, and the finite-branch score must
converge to the dominant eigenspace at the saturated boundary.
"""

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
    left_mode_factors = []
    right_mode_factors = []
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
            left_mode_factors.append(left_factors)
            right_mode_factors.append(right_factors)
            coherence.append(left_factors * right_factors)
        else:
            left_mode_factors.append(np.zeros(len(modes)))
            right_mode_factors.append(np.zeros(len(modes)))
            coherence.append(np.zeros(len(modes)))

    return {
        "left": np.asarray(left_partials),
        "right": np.asarray(right_partials),
        "independent": np.asarray(independent),
        "left_factors": np.asarray(left_mode_factors),
        "right_factors": np.asarray(right_mode_factors),
        "coherence": np.asarray(coherence),
        "invariant_joint": np.asarray(invariant_joint),
    }


def maximum_conditional_mean(probability, values, side):
    """Return max |E[values | one subtree pattern]| under the exact null."""
    conditional_sums = {}
    conditional_masses = {}
    pattern_slice = slice(0, 2) if side == "left" else slice(2, 4)
    for index, pattern in enumerate(PATTERNS):
        key = pattern[pattern_slice]
        conditional_sums.setdefault(key, np.zeros(values.shape[1]))
        conditional_masses.setdefault(key, 0.0)
        conditional_sums[key] += probability[index] * values[index]
        conditional_masses[key] += probability[index]
    conditional_means = [
        conditional_sums[key] / conditional_masses[key]
        for key in conditional_sums
    ]
    return float(np.max(np.abs(conditional_means)))


def check_tangent_geometry(eigenvalues, eigenvectors, inverse, pi, modes):
    """Check normal-space membership and its local information spectrum."""
    category = category_terms(eigenvalues, eigenvectors, inverse, pi, modes, 1.0)
    null_probability = category["independent"]
    null_probability /= null_probability.sum()
    coherence = category["coherence"]

    left_conditional_error = maximum_conditional_mean(
        null_probability, coherence, "left"
    )
    right_conditional_error = maximum_conditional_mean(
        null_probability, coherence, "right"
    )
    conditional_error = max(left_conditional_error, right_conditional_error)
    if conditional_error > 2e-12:
        raise AssertionError(
            "Coherence score is not orthogonal to the independence tangent: "
            f"conditional error={conditional_error}"
        )

    sigma_left = np.einsum(
        "s,si,sj->ij",
        null_probability,
        category["left_factors"],
        category["left_factors"],
    )
    sigma_right = np.einsum(
        "s,si,sj->ij",
        null_probability,
        category["right_factors"],
        category["right_factors"],
    )
    information = np.einsum("s,si,sj->ij", null_probability, coherence, coherence)
    factorization_error = float(
        np.max(np.abs(information - sigma_left * sigma_right))
    )
    if factorization_error > 2e-12:
        raise AssertionError(
            "Tangent information does not equal Sigma_left o Sigma_right: "
            f"error={factorization_error}"
        )

    information_eigenvalues = np.linalg.eigvalsh(information)
    if information_eigenvalues[0] <= 1e-10:
        raise AssertionError(
            "Tangent information is rank deficient: "
            f"eigenvalues={information_eigenvalues}"
        )
    return conditional_error, factorization_error, information_eigenvalues


def check_saturated_tangent(eigenvalues, eigenvectors, inverse, pi, modes):
    """Check that the finite-branch secant converges to the null tangent."""
    category = category_terms(eigenvalues, eigenvectors, inverse, pi, modes, 1.0)
    null_probability = category["independent"]
    null_probability /= null_probability.sum()
    coherence = category["coherence"]

    dominant_eigenvalue = float(np.max(eigenvalues[modes]))
    tolerance = max(1e-8, abs(dominant_eigenvalue) * 1e-6)
    tangent_weights = np.array(
        [
            1.0 if abs(eigenvalues[mode] - dominant_eigenvalue) <= tolerance else 0.0
            for mode in modes
        ]
    )
    faster_gaps = [
        dominant_eigenvalue - eigenvalues[mode]
        for mode, weight in zip(modes, tangent_weights)
        if weight == 0.0
    ]
    branch_length = 32.0 / min(faster_gaps) if faster_gaps else 1.0
    secant_weights = np.exp(
        (eigenvalues[modes] - dominant_eigenvalue) * branch_length
    )
    secant_score = coherence @ secant_weights
    tangent_score = coherence @ tangent_weights
    tangent_error = math.sqrt(
        float(null_probability @ ((secant_score - tangent_score) ** 2))
    )
    if tangent_error > 2e-12:
        raise AssertionError(
            "Finite-branch score did not converge to the saturated tangent: "
            f"L2 error={tangent_error}"
        )
    return tangent_error, tangent_weights


def check_repeated_eigenspace_rotation(
    eigenvalues, eigenvectors, inverse, pi, modes
):
    """Check basis invariance inside a repeated non-stationary eigenspace."""
    repeated_pairs = [
        (first, second)
        for first in range(len(modes))
        for second in range(first + 1, len(modes))
        if abs(eigenvalues[modes[first]] - eigenvalues[modes[second]]) <= 1e-10
    ]
    if not repeated_pairs:
        raise AssertionError("Exact tangent fixture lacks a repeated eigenspace")
    first, second = repeated_pairs[0]
    angle = 0.37
    rotation = np.array(
        [[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]]
    )
    rotated_eigenvectors = eigenvectors.copy()
    rotated_inverse = inverse.copy()
    columns = [modes[first], modes[second]]
    rotated_eigenvectors[:, columns] = eigenvectors[:, columns] @ rotation
    rotated_inverse[columns, :] = rotation.T @ inverse[columns, :]

    original = category_terms(eigenvalues, eigenvectors, inverse, pi, modes, 1.0)
    rotated = category_terms(
        eigenvalues, rotated_eigenvectors, rotated_inverse, pi, modes, 1.0
    )
    branch_length = 0.5
    weights = np.exp(eigenvalues[modes] * branch_length)
    original_score = original["coherence"] @ weights
    rotated_score = rotated["coherence"] @ weights
    rotation_error = float(np.max(np.abs(original_score - rotated_score)))
    if rotation_error > 2e-12:
        raise AssertionError(
            "Weighted score changed under a repeated-eigenspace rotation: "
            f"error={rotation_error}"
        )
    return rotation_error


def information_pseudoinverse(information, relative_cutoff=1e-10):
    """Return a spectrally truncated pseudoinverse and its effective rank."""
    eigenvalues, eigenvectors = np.linalg.eigh(information)
    cutoff = relative_cutoff * float(np.max(eigenvalues))
    retained = eigenvalues > cutoff
    pseudoinverse = (
        eigenvectors[:, retained]
        @ np.diag(1.0 / eigenvalues[retained])
        @ eigenvectors[:, retained].T
    )
    return pseudoinverse, int(np.count_nonzero(retained)), cutoff


def information_decomposition(mean, information, path_weights, pseudoinverse):
    """Split total local information into path and orthogonal components."""
    total = float(mean @ pseudoinverse @ mean)
    path_variance = float(path_weights @ information @ path_weights)
    if path_variance <= 0.0:
        raise AssertionError(f"Non-positive weighted-path variance: {path_variance}")
    path = float((path_weights @ mean) ** 2 / path_variance)
    residual = total - path
    if residual < 0.0 and residual > -2e-12:
        residual = 0.0
    return total, path, residual


def check_omnibus_decomposition(eigenvalues, eigenvectors, inverse, pi, modes):
    """Check the weighted path plus orthogonal-residual enhancement."""
    category = category_terms(eigenvalues, eigenvectors, inverse, pi, modes, 1.0)
    null_probability = category["independent"]
    null_probability /= null_probability.sum()
    coherence = category["coherence"]
    information = np.einsum("s,si,sj->ij", null_probability, coherence, coherence)
    pseudoinverse, effective_rank, cutoff = information_pseudoinverse(information)
    if effective_rank != len(modes):
        raise AssertionError(
            "Exact tangent fixture unexpectedly lost information rank: "
            f"rank={effective_rank}, modes={len(modes)}, cutoff={cutoff}"
        )

    branch_length = 0.5
    persistence = np.exp(eigenvalues[modes] * branch_length)
    on_path_mean = information @ persistence
    on_total, on_path, on_residual = information_decomposition(
        on_path_mean, information, persistence, pseudoinverse
    )
    if abs(on_total - on_path) > 2e-12 or abs(on_residual) > 2e-12:
        raise AssertionError(
            "Correctly specified CTMC signal contains orthogonal residual information: "
            f"total={on_total}, path={on_path}, residual={on_residual}"
        )

    # The zero-length persistence vector is a valid CTMC spectral direction,
    # but it is off the path selected by the fitted length above.  The omnibus
    # statistic should retain information that the fixed directional score
    # cannot use.
    off_path_direction = np.ones(len(modes))
    off_path_mean = information @ off_path_direction
    off_total, off_path, off_residual = information_decomposition(
        off_path_mean, information, persistence, pseudoinverse
    )
    if off_residual <= 1e-6 or off_total <= off_path:
        raise AssertionError(
            "Omnibus residual did not detect an off-path spectral direction: "
            f"total={off_total}, path={off_path}, residual={off_residual}"
        )

    # Exercise the eigenvalue cutoff itself with a controlled nearly singular
    # information matrix.  Exactly one artificial weak direction must be
    # removed from the effective tangent space.
    spectrum, basis = np.linalg.eigh(information)
    spectrum[0] = spectrum[-1] * 1e-14
    near_singular = basis @ np.diag(spectrum) @ basis.T
    truncated, truncated_rank, _ = information_pseudoinverse(near_singular)
    projection = near_singular @ truncated @ near_singular
    projection_error = float(np.max(np.abs(projection - near_singular)))
    if truncated_rank != len(modes) - 1:
        raise AssertionError(
            "Tangent-information cutoff retained the wrong rank: "
            f"rank={truncated_rank}"
        )
    if projection_error > 2e-12:
        raise AssertionError(
            "Truncated tangent pseudoinverse failed its projection identity: "
            f"error={projection_error}"
        )
    return (
        effective_rank,
        on_total,
        on_residual,
        off_total,
        off_path,
        off_residual,
        projection_error,
    )


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

    optimal_weights = np.linalg.solve(covariance, alternative_mean)
    optimality_error = float(np.max(np.abs(optimal_weights - persistence)))
    if optimality_error > 2e-12:
        raise AssertionError(
            "The covariance-adjusted alternative mean did not recover persistence: "
            f"error={optimality_error}"
        )

    def noncentrality(weights):
        return float(
            weights @ alternative_mean
            / math.sqrt(weights @ covariance @ weights)
        )

    persistence_ncp = noncentrality(persistence)
    optimal_ncp = noncentrality(optimal_weights)
    if abs(persistence_ncp - optimal_ncp) > 2e-12:
        raise AssertionError(
            "Persistence weights do not attain the exact linear optimum: "
            f"persistence={persistence_ncp}, optimum={optimal_ncp}"
        )
    return mean_error, optimality_error, persistence_ncp


def main():
    q, pi = build_q(MODEL)
    eigenvalues, eigenvectors, inverse = reversible_eigendecomposition(q, pi)
    zero = int(np.argmin(np.abs(eigenvalues)))
    modes = [index for index in range(len(eigenvalues)) if index != zero]

    conditional_error, factorization_error, information_eigenvalues = (
        check_tangent_geometry(eigenvalues, eigenvectors, inverse, pi, modes)
    )
    tangent_error, tangent_weights = check_saturated_tangent(
        eigenvalues, eigenvectors, inverse, pi, modes
    )
    rotation_error = check_repeated_eigenspace_rotation(
        eigenvalues, eigenvectors, inverse, pi, modes
    )
    (
        effective_rank,
        on_total,
        on_residual,
        off_total,
        off_path,
        off_residual,
        projection_error,
    ) = check_omnibus_decomposition(
        eigenvalues, eigenvectors, inverse, pi, modes
    )
    mean_error, optimality_error, persistence_ncp = check_linear_optimum(
        eigenvalues, eigenvectors, inverse, pi, modes
    )
    print(
        "tangent geometry: "
        f"conditional error={conditional_error:.3g}; "
        f"factorization error={factorization_error:.3g}; "
        f"information eigenvalues={np.array2string(information_eigenvalues, precision=8)}"
    )
    print(
        "saturated tangent: "
        f"L2 error={tangent_error:.3g}; weights={tangent_weights}; "
        f"rotation error={rotation_error:.3g}"
    )
    print(
        "omnibus decomposition: "
        f"rank={effective_rank}; on-path total={on_total:.8g}; "
        f"on-path residual={on_residual:.3g}; off-path total={off_total:.8g}; "
        f"off-path directional={off_path:.8g}; residual={off_residual:.8g}; "
        f"pseudoinverse error={projection_error:.3g}"
    )
    print(
        "linear optimum: "
        f"mean error={mean_error:.3g}; optimality error={optimality_error:.3g}; "
        f"persistence NCP={persistence_ncp:.8g}"
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
