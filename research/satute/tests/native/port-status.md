# SatuTe Phase-1 Port Status

This document records what is implemented and verified in the IQ-TREE native
`--satute` branch. It is deliberately narrower than a manuscript: it is a
code-facing checklist for the Python-reference-first, then native-port workflow.
For a stricter requirement-by-requirement completion audit, see
[completion-audit.md](completion-audit.md).

## Assumptions

- The native analysis runs after IQ-TREE has fitted or accepted the requested
  tree, substitution model and rate-heterogeneity model.
- Phase 1 targets single-alignment reversible models that expose IQ-TREE's
  reversible likelihood kernel and model eigensystem.
- Rate categories are owned by IQ-TREE. Native SatuTe does not infer them from
  scratch.
- The published baseline statistic is the dominant non-stationary eigenspace.
  The promoted `eigenvalue_weighted` statistic is the primary current result;
  `dominant` remains available as the legacy comparator.

## Implemented Native Surface

- `--satute` enables branch-level SatuTe after the main phylogenetic analysis.
- `--satute-alpha NUM` sets the unadjusted significance threshold.
- `--satute-edges FILE` restricts output to branch IDs listed in the file.
- Native output files are:
  - `.sat.stat`: tab-separated branch statistics.
  - `.sat.tree`: Newick tree annotated with the pooled `eigenvalue_weighted`
    statistic and an explicit `satFormula=eigenvalue_weighted` attribute.
  - `.sat.tree.nex`: NEXUS tree for visual inspection.
  - `.sat.branch`: tree with internal branch IDs.

The `.sat.stat` table emits one pooled row per branch and formula. For discrete
rate models it also emits category rows for `dominant` and
`eigenvalue_weighted`; `mixture_likelihood_weighted` is inherently pooled.

## Formula

For each side of a tested branch, ordinary side partial likelihoods are
converted to posterior state probabilities:

```text
P_side(s) = L_side(s) pi_s / sum_x L_side(x) pi_x
```

For a selected set of non-stationary modes, the side coordinate is:

```text
F_side,k = sum_s v_{s,k} P_side(s)
```

The emitted per-site coherence is:

```text
C_i = sum_k w_k F_left,k F_right,k
```

Formula choices:

```text
dominant:             k in the dominant non-stationary eigenspace, w_k = 1
eigenvalue_weighted:  k in all non-stationary modes,
                      w_k = exp((lambda_k - lambda_*) t_c)
mixture_likelihood_weighted:
                      k in all non-stationary modes,
                      w_kc proportional to exp(lambda_k r_c t), with one
                      global log shift and soft null responsibility; a
                      zero-rate invariant component enters only the null
                      denominator
```

Here `lambda_*` is the dominant non-stationary eigenvalue and `t_c` is the tested
branch length after IQ-TREE category-rate rescaling. For homogeneous models,
`t_c` is the branch length. For rate-category rows:

```text
t_c = branch_length * IQTREE_rate_c
```

This means the eigenvectors define the contrast coordinates, while the
eigenvalues determine how much each contrast should still contribute across the
focal branch.

The mixture formula computes the final statistic at site level and estimates
its variance directly across sites. It therefore retains category uncertainty
and between-category variance rather than pooling maximum-posterior category
summaries.

## Rate-Category Contract

Native SatuTe uses IQ-TREE's own empirical-Bayes site-rate pathway:

- `computePatternRates(pattern_rates, pattern_cat)` assigns each pattern to the
  maximum-posterior rate category.
- `getRate(c)` provides the fitted category multiplier.
- The category numbering follows IQ-TREE's `.rate` writer. With `+I+G` or
  `+I+R`, category `0` is the invariant category and positive categories map to
  IQ-TREE's variable-rate classes.

Pure `+I` is intentionally rejected for category-wise SatuTe output because the
IQ-TREE pathway does not emit per-site category rows for that model. `+I`
combined with a discrete rate model is supported when IQ-TREE provides pattern
categories.

## Current Verification

Run all phase-1 checks:

```bash
research/satute/tests/native/run_satute_phase1_checks.sh ./build/iqtree3
```

The bundle verifies:

- JC and GTR informative/saturated smoke cases.
- `--satute-edges` restricts output to the requested branch ID.
- `--satute-edges` rejects absent branch IDs instead of silently analyzing a
  different branch.
- `--satute-alpha` is reflected in `.sat.stat`.
- Rate-category rows match IQ-TREE `.rate` output for `+G`, `+R`, `+I+G` and
  `+I+R`.
- Pure `+I` fails with the expected explicit category-assignment message.
- JC `+G4` category rows match an independent Python numeric reference.
- GTR `+G4` category rows match an independent Python numeric reference for
  `dominant` and `eigenvalue_weighted`; the pooled soft-mixture row also matches
  an independent implementation using the `.iqtree` category table.
- Native pooled `dominant` and `eigenvalue_weighted` rows match the
  independent Python reference for `satC`, `satVar`, `satSE`, `satZ`, `satP`,
  `Decision` and the number of modes.
- The standalone Python reference CLI parses fixed Newick trees, finds the
  requested branch by split, and recomputes homogeneous GTR rows and GTR+G4
  `.rate` category rows from alignment, tree and model inputs. Its
  `--compare-sat-stat` mode fails directly when native `.sat.stat` rows diverge,
  its `--split` option rejects taxon sets that are not a branch in the tree, and
  its `--branch-id` option targets exact native branch IDs through `.sat.stat`.
  The verification bundle includes a negative `--split A,C` check and a
  six-taxon GTR+G4 split plus branch-ID check.

A local verification command used during this audit was:

```bash
research/satute/tests/native/run_satute_phase1_checks.sh ./build/iqtree3 \
  /tmp/iqtree-satute-phase1-checks-missing-edge-20260619
```

It passed together with:

```bash
cmake --build build -j 4
git diff --check
```

## Current Limitations

- Phase 1 rejects supertrees, mixture models and ascertainment-bias correction.
- The native implementation uses scalar branch-side pruning for clarity instead
  of wiring directly into IQ-TREE's vectorized partial-likelihood buffers.
- Protein, codon, morphology, partitioned and Lie-Markov cases have not been
  accepted as verified phase-1 targets.
- The promoted `eigenvalue_weighted` statistic is implemented, independently
  cross-checked and null-calibrated on the current simulation surface. The
  Python comparison retains `all_unweighted` as a non-native research baseline.
- External manuscript/reanalysis files are not synchronized by this checkout
  unless they are made writable in the active workspace.
