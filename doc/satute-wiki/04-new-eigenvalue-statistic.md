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

For GTR, the statistic changes the aggregation step. The published statistic
tests the slowest-decaying direction. The eigenvalue-weighted statistic asks
whether additional model-defined directions still have enough expected
persistence to contribute on the fitted branch. The comparison is therefore:

```text
dominant:             published SatuTe statistic
eigenvalue_weighted:  SatuTe coherence coefficients weighted by relative decay
```

The biological 16S rRNA Tree-of-Life rerun uses this two-method comparison on
the manuscript GTR+F+G4 dataset. In the Eukaryota branch sliding-window
analysis, the curated comparison table reports 895 saturated windows for the
native dominant calculation and 699 saturated windows for the
eigenvalue-weighted statistic out of 1,912 windows. For the yeast branch, the
counts are 83 and 29 saturated windows, respectively. These values are stored
in:

```text
doc/satute-wiki/results/biological/enhanced_sliding_window_rerun_20260619_150104/
```

The independent biological rerun adds both Tree-of-Life examples from the
manuscript. The protein analysis uses `protein_based_2D_tree` and the 16S rRNA
analysis uses `rRNA_based_3D_tree`. In the protein 2D tree, the published
SatuTe rerun marks 138 and 54 saturated windows for the Eukaryota and yeast
branches; the eigenvalue-weighted calculation marks 5 and 8. In the rRNA 3D
tree, the corresponding counts are 822 and 70 for published SatuTe, and 699
and 29 for the eigenvalue-weighted calculation. These values are stored in:

```text
doc/satute-wiki/results/biological/independent_enhanced_sliding_window_20260619_151135/
```

The simulation grid is still required for calibration. The biological reruns
are empirical examples showing how the statistic changes real sliding-window
analyses; they do not replace type I error and power evaluation under the
original SatuTe simulation design.
