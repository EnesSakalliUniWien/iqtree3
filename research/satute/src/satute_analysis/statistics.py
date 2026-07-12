"""Small dependency-free statistical helpers for paired experiments."""

import math

import numpy as np


def wilson_interval(successes, total, z=1.96):
    fraction = successes / total
    denominator = 1.0 + z * z / total
    center = (fraction + z * z / (2.0 * total)) / denominator
    half = z * math.sqrt(
        fraction * (1.0 - fraction) / total + z * z / (4.0 * total * total)
    ) / denominator
    return center - half, center + half


def mean_interval(values, z=1.96):
    values = np.asarray(values, dtype=float)
    mean = float(np.mean(values))
    half = z * float(np.std(values, ddof=1)) / math.sqrt(len(values))
    return mean, mean - half, mean + half


def paired_bootstrap_interval(refined, legacy, rng, replicates=4000):
    differences = np.asarray(refined, dtype=int) - np.asarray(legacy, dtype=int)
    counts = np.array(
        [np.sum(differences == -1), np.sum(differences == 0), np.sum(differences == 1)]
    )
    draws = rng.multinomial(len(differences), counts / len(differences), size=replicates)
    bootstrap = (draws[:, 2] - draws[:, 0]) / len(differences)
    return float(np.quantile(bootstrap, 0.025)), float(np.quantile(bootstrap, 0.975))
