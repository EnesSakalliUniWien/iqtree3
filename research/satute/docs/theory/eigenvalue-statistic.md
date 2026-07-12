# Eigenvalue-Weighted Statistic

The published SatuTe statistic uses the dominant non-stationary eigenspace.
Under JC, the three non-stationary eigenvalues are equal, so the dominant
eigenspace contains all three non-stationary directions. Under a typical GTR
model, the non-stationary eigenvalues differ, and the dominant eigenspace often
contains only the slowest-decaying direction.

The eigenvalue-weighted version keeps the SatuTe coherence coefficients but
uses the fitted rate matrix to weight non-stationary directions by their
expected persistence across the focal branch:

```text
C_site_lambda = sum_{k != stationary}
                exp((lambda_k - lambda_*) * branch_length)
                * F_left(k) * F_right(k)
```

Here `lambda_*` is the dominant non-stationary eigenvalue. The relative weight
is one for the dominant mode and smaller for faster-decaying modes. The
eigenvectors define the directions being compared; the eigenvalues define the
expected persistence of those directions across the focal branch.

For JC, all non-stationary eigenvalues are equal. The weights are therefore one
for every non-stationary direction, and the weighted statistic gives the same
normalized test as the published SatuTe statistic.

## Null Hypothesis

The weighted statistic keeps the SatuTe null hypothesis: the two subtree
patterns on the sides of the focal branch are independent. Write the mode
factors for one site as:

```text
X_i = < L_left / P_left, h_i >
Y_i = < L_right / P_right, h_i >
```

Then the ordinary SatuTe mode coefficient is:

```text
C_i = X_i * Y_i
```

and the weighted coefficient is the deterministic linear combination:

```text
C_lambda(t) = sum_i w_i(t) * X_i * Y_i
```

Under subtree independence, the left and right factors are independent. Since
each non-stationary eigenvector is orthogonal to the stationary component, the
same argument used by SatuTe gives:

```text
E[X_i] = 0
E[Y_i] = 0
```

Therefore:

```text
E[C_i] = E[X_i * Y_i] = E[X_i] * E[Y_i] = 0
E[C_lambda(t)] = sum_i w_i(t) * E[C_i] = 0
```

So weighting does not move the null mean. The null variance is obtained by
expanding the weighted square:

```text
Var(C_lambda(t))
  = sum_i sum_j w_i(t) * w_j(t)
      * sigma_A(i,j) * sigma_B(i,j)
```

where:

```text
sigma_A(i,j) = E[X_i * X_j]
sigma_B(i,j) = E[Y_i * Y_j]
```

This is exactly the SatuTe variance structure with the mode weights included.
For independent sites, the sample mean has variance `Var(C_lambda(t)) / n`.
With the plug-in variance estimator, the central limit theorem gives the same
one-sided normal null test used by SatuTe:

```text
Z_lambda(t) = sqrt(n) * mean(C_lambda(t)) / sigma_lambda(t)
```

Thus the weighted method is mathematically a SatuTe-compatible linear spectral
test. It changes which non-stationary directions enter the statistic, but it
does not change the null hypothesis or the null calibration argument. As in
SatuTe, maximum-likelihood topology selection is a separate source of
double-dipping, so the Bonferroni correction is retained for that scenario.

## Null Calibration Smoke Check

The corrected INDELible development check uses the same GTR matrix for
simulation and evaluation. The translation matters because IQ-TREE writes DNA
frequencies as `A,C,G,T` and GTR rates as `AC,AG,AT,CG,CT,GT`, while INDELible
expects frequencies as `T,C,A,G` and GTR rates as `CT,AT,GT,AC,CG,AG`.

A near-null smoke check simulated 500 alignments with 1,000 sites, the
five-taxon tree, the focal branch length set to 100 substitutions per site, and
the skewed GTR model `GTR_SKEW_BOTH`. The true-tree fixed-length scenario is
the cleanest calibration check because it avoids topology-selection effects.
At nominal `alpha = 0.05`, both the dominant and eigenvalue-weighted formulas
called 26 of 500 replicates informative:

```text
formula               informative / evaluated   fraction   Wilson 95% CI
dominant              26 / 500                   0.052      0.036-0.075
eigenvalue_weighted   26 / 500                   0.052      0.036-0.075
```

The mean null `Z` score was 0.039 and the standard deviation was 0.992 for both
formulas. This supports the expected null centering and normal scale in the
fixed-tree setting. It is still a smoke check: the full calibration must cover
the complete tree-case, branch-length, sequence-length and topology-selection
grid.

For GTR, the statistic changes the aggregation step. The published statistic
tests the slowest-decaying direction. The eigenvalue-weighted statistic asks
whether additional model-defined directions still have enough expected
persistence to contribute on the fitted branch. The comparison is therefore:

```text
published_satute:     original workflow output from the SatuTe analysis
dominant:             native recomputation of the published dominant-mode statistic
eigenvalue_weighted:  native recomputation with coefficients weighted by relative decay
```

The `dominant` and `published_satute` quantities are meant to represent the
same SatuTe formula, but they are not the same provenance. Published SatuTe is
the original workflow output. The dominant calculation is recomputed through
the native IQ-TREE input path. Differences between those two traces are
therefore implementation/input differences; formula effects should be read from
the native `dominant` versus native `eigenvalue_weighted` comparison.

The biological 16S rRNA Tree-of-Life rerun uses this native formula comparison on
the manuscript GTR+F+G4 dataset. In the Eukaryota branch sliding-window
analysis, the curated comparison table reports 895 saturated windows for the
native dominant-mode recomputation and 699 saturated windows for the
eigenvalue-weighted statistic out of 1,912 windows. For the yeast branch, the
counts are 83 and 29 saturated windows, respectively. These values are stored
in:

```text
research/satute/artifacts/releases/biological/enhanced_sliding_window_rerun_20260619_150104/
```

The independent biological rerun adds both Tree-of-Life examples from the
manuscript. The protein analysis uses `protein_based_2D_tree` and the 16S rRNA
analysis uses `rRNA_based_3D_tree`. In the protein 2D tree, the published
SatuTe rerun marks 138 and 54 saturated windows for the Eukaryota and yeast
branches; the eigenvalue-weighted calculation marks 5 and 8. In the rRNA 3D
tree, the corresponding counts are 822 and 70 for published SatuTe, and 699
and 29 for the eigenvalue-weighted calculation. These values are stored in:

```text
research/satute/artifacts/releases/biological/independent_enhanced_sliding_window_20260619_151135/
```

The simulation grid is still required for calibration. The biological reruns
are empirical examples showing how the statistic changes real sliding-window
analyses; they do not replace type I error and power evaluation under the
original SatuTe simulation design.
