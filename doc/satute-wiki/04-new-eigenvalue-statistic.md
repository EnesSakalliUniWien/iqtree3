# Eigenvalue-Weighted Statistic

The published SatuTe statistic uses the dominant non-stationary eigenspace.
Under JC, the three non-stationary eigenvalues are equal, so the dominant
eigenspace contains all three non-stationary directions. Under a typical GTR
model, the non-stationary eigenvalues differ, and the dominant eigenspace often
contains only the slowest-decaying direction.

The unweighted all-mode statistic is:

```text
C_site_all = sum_{k != stationary} F_left(k) * F_right(k)
```

This statistic includes every non-stationary eigenvector. It can detect
agreement in directions not used by the dominant statistic, but it gives equal
weight to modes with different decay rates.

The eigenvalue-weighted version is:

```text
C_site_lambda = sum_{k != stationary}
                exp((lambda_k - lambda_*) * branch_length)
                * F_left(k) * F_right(k)
```

Here `lambda_*` is the dominant non-stationary eigenvalue. The relative weight
is one for the dominant mode and smaller for faster-decaying modes. The
eigenvectors define the directions being compared; the eigenvalues define the
expected persistence of those directions across the focal branch.

In the GTR formula comparison, the `gtr_5.00` case produced different decisions:

```text
dominant:             p = 0.321608  saturated
all_unweighted:       p = 0.011107  informative
eigenvalue_weighted:  p = 0.189025  saturated
```

The additional signal in `all_unweighted` came from a mode with eigenvalue
`-1.53061`. Across the fitted branch, its relative decay weight was `0.001962`.
The slow dominant mode had eigenvalue `-1.02041` and relative decay weight
`0.015672`. The weighted statistic therefore reduced the contribution of the
faster mode instead of treating it as equally persistent.

The current interpretation is limited to method development. `all_unweighted`
is a sensitivity diagnostic for non-dominant spectral signal.
`eigenvalue_weighted` is the model-aware extension to test in the full
simulation grid. Its variance estimator and p-value calibration still need to
be evaluated under the original SatuTe simulation design and under controlled
GTR alternatives.
