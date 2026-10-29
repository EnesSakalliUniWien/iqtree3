#!/usr/bin/env python3
"""Check the sliding-window driver against the published rolling construction.

The reference here is a direct transcription of the pandas/numpy code that
produced the curated biological releases: a rolling weighted mean over the
coherence series, divided by sqrt(variance / window), then shifted back by
half a window. The driver computes the same thing without pandas, so this
pins the port rather than the formula.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "experiments" / "005_biological_windows"))

from sliding_window import window_scores, window_weights  # noqa: E402


def reference_scores(coherence, variance, window_size):
    weights = np.linspace(1, 0.5, window_size)

    def score(values):
        return np.average(values, weights=weights) / np.sqrt(variance / window_size)

    rolled = (
        pd.Series(coherence)
        .rolling(window=window_size)
        .apply(score, raw=True)
        .shift(-(window_size // 2))
    )
    return list(rolled)


def check_matches_reference():
    rng = np.random.default_rng(20261029)
    for window_size in (4, 7, 36):
        for length in (window_size, window_size + 1, 120):
            coherence = list(rng.normal(0.3, 1.1, size=length))
            variance = 0.87
            mine = window_scores(coherence, variance, window_size)
            theirs = reference_scores(coherence, variance, window_size)
            assert len(mine) == len(theirs) == length
            for index, (got, want) in enumerate(zip(mine, theirs)):
                if got is None:
                    assert np.isnan(want), f"w={window_size} n={length} i={index}"
                else:
                    assert abs(got - want) < 1e-12, (
                        f"w={window_size} n={length} i={index}: {got} vs {want}"
                    )


def check_padding_matches_release():
    # The frozen protein release holds 2,596 rows for a 36-site window: 17
    # blank, 2,561 scores, 18 blank. Reproduce that shape exactly.
    scores = window_scores([0.5] * 2596, 1.0, 36)
    leading = next(i for i, value in enumerate(scores) if value is not None)
    trailing = next(i for i, value in enumerate(reversed(scores)) if value is not None)
    finite = sum(1 for value in scores if value is not None)
    assert (leading, finite, trailing) == (17, 2561, 18), (leading, finite, trailing)


def check_weights():
    weights = window_weights(36)
    assert len(weights) == 36
    assert abs(weights[0] - 1.0) < 1e-12
    assert abs(weights[-1] - 0.5) < 1e-12
    assert all(a > b for a, b in zip(weights, weights[1:]))


def check_rejects_bad_variance():
    for bad in (0.0, -1.0):
        try:
            window_scores([1.0] * 40, bad, 36)
        except ValueError:
            continue
        raise AssertionError(f"accepted satVar={bad}")


def main():
    check_weights()
    check_matches_reference()
    check_padding_matches_release()
    check_rejects_bad_variance()
    print("Sliding-window driver: passed")


if __name__ == "__main__":
    main()
