# SatuTe Phase-1 Completion Audit

This audit maps the active implementation objective to evidence in the current
checkout. It is intentionally stricter than a status note: a requirement is
marked complete only when the repository contains the implementation and an
explicit verifier covers the behavior.

## Objective Interpreted For This Checkout

The working objective is:

```text
Python reanalysis tool first
Native IQ-TREE --satute first
Python reference first, then port to IQ-TREE
```

The interpretation used by this checkout is:

1. Build an independent Python reference/reanalysis path that can recompute the
   SatuTe formula outside native IQ-TREE for controlled validation cases.
2. Port the same formula into IQ-TREE as native `--satute` output.
3. Verify the native rows against the Python reference, including GTR and
   IQ-TREE-owned rate categories.

This now means a general fixed-tree nucleotide Python reference for the phase-1
model surface, not a replacement for every upstream SatuTe workflow.

## Requirement Audit

| Requirement | Current evidence | Status |
| --- | --- | --- |
| Independent Python reference exists before trusting native rows | `satute_reference.py` parses fixed Newick trees, finds a requested split or native branch ID, recomputes homogeneous JC/GTR and IQ-TREE `.rate` category rows; `compare_formulas.py` independently recomputes native pooled rows. | Complete for nucleotide fixed-tree phase-1 validation. |
| Native IQ-TREE command-line surface exists | `--satute`, `--satute-alpha`, and `--satute-edges` are parsed in `utils/tools.cpp`; `PhyloTree::computeSatuTe` is called after the main analysis in `main/phyloanalysis.cpp`. | Complete for phase 1. |
| Native implementation computes the posterior/eigenmode statistic | `tree/satute.cpp` converts side partials to posterior state probabilities, projects them onto model right eigenvectors, computes per-site coherence, and emits `satC`, variance, standard error, z-score and p-value. | Complete for reversible stationary models that pass the current gates. |
| Formula variants are emitted | Native output includes `dominant`, `eigenvalue_weighted`, and the pooled `mixture_likelihood_weighted`. The latter uses soft null responsibilities, globally rescaled `exp(lambda_k r_c t)` persistence weights, the correct invariant joint likelihood for `+I`, and direct site-level variance. The Python comparison keeps `all_unweighted` as an offline baseline. | Complete. |
| JC multiplicity behaves correctly | The formula comparison verifies that JC emits identical pooled rows for `dominant` and `eigenvalue_weighted`. | Complete for tested JC cases. |
| GTR directionality is represented | GTR examples use a non-uniform stationary distribution and asymmetric exchangeability choices; the dominant row can use a single dominant non-stationary mode while `eigenvalue_weighted` compares all non-stationary modes. | Complete for tested GTR cases. |
| Rate categories come from IQ-TREE, not a custom SatuTe estimator | Native SatuTe calls `site_rate->computePatternRates(pattern_rates, pattern_cat)` and `site_rate->getRate(c)`. The regression compares SatuTe category counts and rates to IQ-TREE `.rate` output. | Complete for `+G`, `+R`, `+I+G`, and `+I+R` cases covered by the test. |
| Pure `+I` is handled honestly | The rate-category regression requires pure `+I` to fail with an explicit message because IQ-TREE does not emit per-site category rows for that model. | Complete. |
| Branch subset behavior is checked | The smoke test verifies `--satute-edges`, `--satute-alpha`, and failure on an absent branch ID. | Complete. |
| Native rows are checked against Python | The phase-1 bundle compares native rows to independent Python rows for `satC`, `satVar`, `satSE`, `satZ`, `satP`, decision and mode count. It includes four-taxon checks and a six-taxon GTR+G4 split plus native branch-ID targeting. | Complete for the covered nucleotide phase-1 models. |
| Build and whitespace hygiene pass | Current audit command set includes `research/satute/tests/native/run_satute_phase1_checks.sh`, `cmake --build build -j 4`, and `git diff --check`. | Complete for the current local checkout. |

## Current Verification Command

The current local verification command is:

```bash
research/satute/tests/native/run_satute_phase1_checks.sh ./build/iqtree3 \
  /tmp/iqtree-satute-phase1-checks
cmake --build build -j 4
git diff --check
```

That command set passed in the current checkout.

## Remaining Phase-1 Boundaries

The phase-1 IQ-TREE port and Python reference are verified for the stated
nucleotide fixed-tree surface. These boundaries remain explicit:

- Protein, codon, morphology, partitioned, mixture, Lie-Markov and
  ascertainment-corrected analyses are not accepted phase-1 targets.
- The promoted `eigenvalue_weighted` statistic is numerically implemented,
  independently cross-checked and covered by exact-null and paired simulation
  calibration; `dominant` remains the published legacy comparator.
- External manuscript and reanalysis artifacts outside this checkout are outside
  this repository audit.

## Next Completion Step

The next completion step is packaging: make the biological sliding-window
analysis self-contained in this repository and preserve a publication rerun
manifest for the promoted `eigenvalue_weighted` statistic.
