# SatuTe Formula Comparison

This note compares the two native posterior/eigenvector formulas plus one
offline all-mode baseline on simulated four-taxon datasets.

## Formulas

### Dominant-Posterior Formula

This is the first formula emitted by `tree/satute.cpp` and documented in
`README.md`. It is also the formula used for the annotations in `.sat.tree`.

For each side of a tested branch, ordinary partial likelihoods are converted
to posterior state probabilities:

```text
P_side(s) = L_side(s) * pi(s) / sum_x L_side(x) * pi(x)
```

For the dominant non-zero eigenmode or degenerate dominant eigenspace:

```text
F_side(k) = sum_s right_eigenvector(s,k) * P_side(s)
C_site    = sum_k F_left(k) * F_right(k)
```

This is the SatuTe-compatible statistic: it focuses on the slowest-decaying
non-stationary mode, which is the mode that controls long-time branch
saturation.

### All-Mode Posterior Formula

This is an offline unweighted eigenvector comparator used by the independent
Python analysis; it is not emitted by native `--satute`:

```text
C_site_all = sum_{k != stationary} F_left(k) * F_right(k)
```

It uses the same posterior construction, but includes every non-zero
eigenmode, not only the dominant non-zero eigenspace.

### Eigenvalue-Decay Weighted Formula

This is the eigenvector-plus-eigenvalue comparator now emitted by native
`--satute`:

```text
C_site_decay = sum_{k != stationary}
               exp((lambda_k - lambda_1) * branch_length)
               * F_left(k) * F_right(k)
```

The eigenvectors define the signal coordinates. The eigenvalues define how
fast each coordinate decays over the focal branch. More negative eigenvalues
receive smaller relative weights on long branches, while the dominant
non-stationary mode keeps weight one.

For rate-category models, `branch_length` in this expression is the tested
branch length multiplied by the IQ-TREE category rate. Native `--satute` obtains
the maximum-posterior pattern category and category multiplier from IQ-TREE,
then emits category-specific weights and a pooled row whose variance uses the
SatuTe global-variance rule. Pure `+I` is excluded from category-wise output
unless IQ-TREE combines it with a discrete rate model that provides
`computePatternRates` categories.

## Reproducible Command

Set `PYTHON_BIN` if the default `python3` does not have NumPy:

```bash
"${PYTHON_BIN:-python3}" research/satute/tests/native/compare_formulas.py \
  --iqtree ./build/iqtree3 \
  --outdir /tmp/iqtree-satute-formula-compare
```

## Results

Simulation used 5,000 sites on the split `AB|CD`, with terminal branch lengths
`0.05`. GTR examples used:

```text
GTR{1,2,1,1,2,1}+F{0.30,0.20,0.20,0.30}
```

| Case | Model | Simulated branch | Optimized branch | Dominant p/decision | All-mode p/decision | Decay-weighted p/decision |
|---|---|---:|---:|---|---|---|
| `jc_0.50` | JC | 0.50 | 0.50871 | 0 / informative | 0 / informative | 0 / informative |
| `jc_8.00` | JC | 8.00 | 10.00000 | 0.612926 / saturated | 0.612926 / saturated | 0.612926 / saturated |
| `gtr_0.50` | GTR | 0.50 | 0.497267 | 0 / informative | 0 / informative | 0 / informative |
| `gtr_4.00` | GTR | 4.00 | 3.19089 | 0.00896166 / informative | 0.00644015 / informative | 0.00416481 / informative |
| `gtr_5.00` | GTR | 5.00 | 4.07273 | 0.321608 / saturated | 0.011107 / informative | 0.189025 / saturated |
| `gtr_8.00` | GTR | 8.00 | 10.00000 | 0.765772 / saturated | 0.962115 / saturated | 0.770021 / saturated |

The independent Python recomputation checks IQ-TREE's native `--satute` output
for the `dominant` and `eigenvalue_weighted` pooled rows. The
`all_unweighted` values in this table are computed only by the independent
comparison code.

## Interpretation

The unweighted all-mode formula reacts more strongly near the GTR transition
region. In the `gtr_5.00` case, it reports significant signal, while both the
`dominant` formula and the `eigenvalue_weighted` formula report
saturation.

The mode decomposition explains why. At `gtr_5.00`, the fitted branch length is
`4.07273`. The two faster modes have eigenvalue `-1.53061` and branch decay
weight `0.001962`; the slow dominant mode has eigenvalue `-1.02041` and branch
decay weight `0.015672`. The unweighted all-mode formula gives the fast modes
equal influence, so one fast mode remains individually significant. The
`eigenvalue_weighted` formula recognizes that this mode should barely survive across
the branch.

This is why eigenvalues are important. Eigenvectors identify what pattern of
state contrast is being measured, but eigenvalues say whether that contrast is
slow enough to matter over the tested evolutionary distance.

For the IQ-TREE integration, the `dominant` formula remains the
published SatuTe baseline because it:

- matches the SatuTe-compatible posterior/eigenvector construction
- focuses on the slowest non-stationary mode that controls saturation
- has the cleanest existing calibration

The decay-weighted version is a native research diagnostic. The offline
unweighted all-mode version shows how additional eigenvector coordinates
can recover signal; the decay-weighted version tests whether that signal is
plausible across the focal branch once the eigenvalues are respected. Before
either replaces the published statistic, its null variance and p-value
calibration need broader simulation validation.
