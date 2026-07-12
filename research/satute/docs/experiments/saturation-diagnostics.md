# Saturation Diagnostics: Implementation Audit and Analysis Contract

## Purpose

The native SatuTe output contains three quantities that answer different
questions:

$$
\bar C(t),\qquad \operatorname{SE}\{\bar C(t)\},\qquad
Z(t)=\frac{\bar C(t)}{\operatorname{SE}\{\bar C(t)\}}.
$$

They must be inspected separately. The raw coherence is a formula-dependent
effect scale, the standard error describes uncertainty, and the standardized
statistic is a detection statistic. None is, by itself, a guaranteed monotonic
measure of evolutionary saturation.

The focused diagnostic pipeline is
[`saturation_diagnostics.py`](saturation_diagnostics.py). It reads the three
native formulas directly from `.sat.stat`, reads IQ-TREE's fitted rate-category
table from `.iqtree`, estimates formula-specific null thresholds on one half of
the saturated simulations, and evaluates those thresholds on the held-out
half.

## Why the current implementation has its present form

### Dominant eigenspace

The original statistic keeps the slowest nonstationary eigenspace. This is a
low-dimensional test of the contrast expected to persist longest. It is also
the native tree-annotation formula, preserving the original SatuTe interface.

### Relative eigenvalue weights

The current `eigenvalue_weighted` formula uses

$$
w_{kc}(t)=\exp\{(\lambda_k-\lambda_*)r_ct\}.
$$

Subtracting the dominant eigenvalue makes the largest weight exactly one. This
prevents a common exponentially small factor from causing numerical underflow
and leaves a category-specific $Z$-score unchanged when both its mean and
standard error receive the same scale factor. It is not innocuous after rate
categories are pooled, because the removed common factor differs by category.

### Hard rate-category assignment

The original and relative-weighted rows use IQ-TREE's
`computePatternRates` maximum-posterior category. This matches the original
workflow, which selected the largest `.siteprob` posterior, and provides clear
category-specific diagnostic rows. Its limitation is that category uncertainty
is discarded.

### Factorized null variance

Under the saturated null, the two sides of the branch are independent. For
mode coordinates $F^{(L)}_k$ and $F^{(R)}_k$, this motivates

$$
\operatorname{E}_0\!\left[
F^{(L)}_kF^{(L)}_jF^{(R)}_kF^{(R)}_j
\right]
=
\operatorname{E}_0\!\left[F^{(L)}_kF^{(L)}_j\right]
\operatorname{E}_0\!\left[F^{(R)}_kF^{(R)}_j\right].
$$

The native variance is the corresponding weighted sum. Category variances are
pooled under the null assumption that the conditional category means are zero.
Finite-sample category selection and parameter estimation can violate the
normal approximation even when this factorization is correct asymptotically.

### Soft-mixture statistic

For positive-rate categories, `mixture_likelihood_weighted` uses the
focal-branch independence likelihoods

$$
b_{sc}=q_cA_{sc}B_{sc}
$$

and soft responsibilities $\tau_{sc}=b_{sc}/\sum_db_{sd}$. Its site score is

$$
u_s(t)=\sum_c\tau_{sc}\sum_{k\ne0}
e^{\lambda_kr_ct}C_{sc,k}.
$$

This satisfies the exact identity

$$
\frac{P_t(\partial_s)}{P_\infty(\partial_s)}=1+u_s(t)
$$

when every category rate is positive. The native calculation applies one
common log-weight shift to all positive-rate categories and modes to prevent
underflow. This multiplies every site score by the same constant and therefore
does not change `satZ`.

With an explicit zero-rate invariant category, the model's infinite-branch
null is not category-wise independence because the invariant transition matrix
remains the identity. The null denominator therefore uses

$$
q_0J_{s0}+\sum_{c:r_c>0}q_cA_{sc}B_{sc},
\qquad
J_{s0}=\sum_a\pi_aL_{A,0}(a)L_{B,0}(a).
$$

Only positive-rate categories enter the numerator. The invariant component is
present in both the finite- and infinite-branch models and hence contributes
zero evidence about branch length. Because the final score is constructed at
site level after integrating categories, its variance is estimated directly
across final site scores.

## Monotonic model-based saturation index

For the diagnostic figure, set $a_k=1$ for every nonstationary eigenmode of the
reversible symmetrized rate matrix and define

$$
I(t)=\sum_cq_c\sum_{k\ne0}\exp(2\lambda_kr_ct),
$$

$$
S(t)=1-\frac{I(t)}{I(0)}.
$$

This choice measures the normalized spectral energy remaining outside the
stationary component. Since $q_c>0$, $r_c\ge0$, and $\lambda_k<0$,

$$
\frac{dI(t)}{dt}
=
\sum_cq_c\sum_{k\ne0}
2\lambda_kr_c\exp(2\lambda_kr_ct)
\le0.
$$

Therefore $S(t)$ is monotonic nondecreasing under a fixed model. It is a model-
based ruler, not a significance test. With an explicit zero-rate invariant
category, its limiting value can remain below one; that behavior should be
reported rather than silently renormalized away.

## Native output contract

The scale is now computed directly by `tree/satute.cpp`. Every `.sat.stat` row
contains:

```text
InformationFraction = I(t) / I(0)
SaturationIndex      = 1 - InformationFraction
```

For category rows, $t$ is replaced by the category effective length $r_ct$.
For pooled rows, spectral energy is averaged using IQ-TREE's full-precision
category proportions:

$$
\frac{I(t)}{I(0)}
=
\frac{1}{\sum_cq_c}
\sum_cq_c\frac{1}{K}
\sum_{k\ne0}\exp(2\lambda_kr_ct),
$$

where $K$ is the number of nonstationary modes. The value is identical across
the three SatuTe formula rows because it depends only on the fitted model,
branch length and rate distribution. NEXUS tree annotations expose the pooled
values as `satInfo` and `satIndex`.

## Figure contract

The generated six-panel figure separates:

1. mean raw coherence with a 95% confidence interval;
2. mean standard error on a logarithmic scale;
3. median standardized statistic with its 10th--90th percentile band;
4. category-specific dominant persistence;
5. calibrated informative probability with Wilson intervals;
6. the monotonic model-based saturation index.

All panels use the same logarithmic branch-length axis. Formula identity uses
color, line style, and marker shape, so it remains readable without color.

## Stepwise assumptions tested

The pipeline writes `assumption_checks.md` and tests:

1. `satZ = satC / satSE` numerically;
2. positive category proportions summing to one and nonnegative rates;
3. exact dominant/relative weight normalization;
4. monotonic category persistence;
5. monotonicity and boundary behavior of $S(t)$;
6. reduction of mean absolute raw coherence from the shortest branch to the
   saturated branch;
7. the standard-normal null assumption using the null mean and standard
   deviation of `satZ`;
8. held-out type-I error after empirical calibration.

The fixed-tree, fixed-model experiment isolates formula behavior. A later
end-to-end experiment must repeat calibration after estimating the model, rate
distribution, topology and branch lengths.

## Verified fixed-model results

The implemented experiment used 1,000-site alignments, 100 replicates at each
of 14 finite branch lengths, and 500 saturated-branch replicates. The latter
were divided equally between threshold training and held-out validation.

The algebraic identities, category contract, relative-weight normalization,
hard-category pooling identities, persistence monotonicity and saturation-index
monotonicity all passed. The maximum discrepancy in
`satZ = satC / satSE` was $2.02\times10^{-8}$; the maximum recomputation errors
for pooled coherence and variance were $8.72\times10^{-10}$ and
$8.48\times10^{-10}$.

The native saturation scale also passed its independent check. Values were
identical across formulas; recomputation from the rounded `.iqtree` rate table
agreed within $1.22\times10^{-5}$. Under homogeneous JC, the native result
matched the closed form $S(t)=1-\exp(-8t/3)$ within
$5.9\times10^{-11}$.

The refined weighted likelihood is checked independently by
[`verify_weighted_likelihood.py`](verify_weighted_likelihood.py). Exact
enumeration recovered the positive-rate mixture likelihood-ratio identity with
maximum error $5.4\times10^{-14}$ and the corrected `+I` identity with maximum
error $1.1\times10^{-11}$ at branch length $10{,}000$. The corresponding null
means were below $2\times10^{-16}$. The native implementation also matched the
independent Python reference for both `GTR+G4` and `JC+I+G4` rows within
$2\times10^{-4}$.

### Simulation-independence correction

An initial version launched each replicate as a separately seeded ALISIM
process. Those short SPRNG streams were correlated and produced a false
positive centering shift. Across 100,000 sites from those processes, the
empirical 256-pattern distribution gave

$$
\chi^2=446.99\qquad(df=255),
$$

and overrepresented high-coherence patterns. This was not a formula failure:
the exact pattern probabilities reproduced IQ-TREE's log likelihood to
$1.2\times10^{-5}$, and a single one-million-site ALISIM stream gave
$\chi^2=235.89$. ALISIM's native `--num-alignments` mode gave
$\chi^2=253.45$ across 100 independent alignments.

The diagnostic runner now generates every replicate set with one independent
`--num-alignments` batch stream. All results below supersede the earlier
separate-process simulation results.

### Standard-normal null verification

With the corrected simulation design, the standard-normal null assumption
passed:

| Formula | Null mean $Z$ | Null SD | Nominal rejection at 1.645 | Held-out calibrated rejection |
|---|---:|---:|---:|---:|
| `dominant` | 0.024 | 0.967 | 0.044 | 0.052 |
| `eigenvalue_weighted` | 0.024 | 0.967 | 0.044 | 0.052 |
| `mixture_likelihood_weighted` | 0.029 | 0.957 | 0.056 | 0.068 |

Exact enumeration of all $4^4=256$ patterns gave null expectations of zero for
the hard statistic and $-3.5\times10^{-18}$ for the normalized soft statistic.
Direct multinomial sampling with 20,000 replicates at each alignment length
gave null rejection rates between 0.0483 and 0.0518.
The exact enumeration and multinomial assertions are reproduced by
[`verify_saturation_null.py`](verify_saturation_null.py).

The native alignment-length experiment used 1,000 batch-generated replicates
per value:

| Sites | Dominant/weighted mean $Z$ | Dominant/weighted SD | Soft mean $Z$ | Soft SD |
|---:|---:|---:|---:|---:|
| 100 | 0.019 | 0.977 | 0.011 | 0.976 |
| 500 | -0.030 | 0.984 | -0.017 | 0.992 |
| 1,000 | -0.007 | 1.007 | -0.012 | 1.009 |
| 5,000 | 0.027 | 1.001 | 0.036 | 1.012 |
| 10,000 | 0.048 | 0.993 | 0.056 | 0.999 |

Every mean confidence interval included zero, and the empirical variance of
the replicate means agreed with the internal variance estimates. No revised
null statistic or empirical centering correction is indicated.

Calibrated informative fractions in the transition region were:

| Branch length | Dominant | Relative weighted | Soft mixture |
|---:|---:|---:|---:|
| 30 | 0.94 | 1.00 | 0.99 |
| 50 | 0.45 | 0.62 | 0.63 |
| 75 | 0.21 | 0.21 | 0.18 |
| 100 | 0.09 | 0.09 | 0.07 |

These values establish behavior only for the fixed GTR+G4 model and tree used
here. At branch length 50, the paired soft-minus-relative difference was 0.01
with a bootstrap 95% interval of $[-0.04,0.06]$, providing no evidence of a
power difference. These results are not a general ranking across phylogenetic
conditions.

### Dedicated `+I+G4` and `+I+R4` paired grid

The invariant-category correction was evaluated separately with
[`weighted_invariant_simulations.py`](weighted_invariant_simulations.py). The
script enumerates all $4^4=256$ site patterns under fixed `GTR+I+G4` and
`GTR+I+R4` models, then draws paired multinomial alignments. Every alignment is
evaluated by both the former counterfactual construction and the refined
construction. The grid contains

$$
n\in\{100,500,1000,5000,10000\},
$$

$$
t\in\{0.5,2,5,10,20,50,200,10000\},
$$

with 6,000 null and 6,000 alternative replicates per cell. Half of each null
sample trains a formula-specific 95th-percentile threshold; the other half is
held out. In total, 960,000 simulated alignments were evaluated by both
formulas.

The exact checks passed for both models and all eight branch lengths. The
refined null mean was below $1.3\times10^{-16}$, and the maximum error in

$$
\frac{p_t(s)}{p_\infty(s)}
=1+e^\delta u_s(t)
$$

was $1.6\times10^{-11}$. In contrast, the legacy counterfactual score retained
a positive expectation under the correct null, approximately $0.0355$ at
$t=10{,}000$ for both mixtures.

This structural shift dominated nominal false-positive control. Across all 80
model-by-branch-by-length cells, the refined rejection rate at the ordinary
$N(0,1)$ cutoff averaged 0.0488 and ranged from 0.0377 to 0.0577. The legacy
rate averaged 0.901 and ranged from 0.500 to 1.000. At $t=10{,}000$, the
comparison was:

| Sites | Legacy `+I+G4` | Refined `+I+G4` | Legacy `+I+R4` | Refined `+I+R4` |
|---:|---:|---:|---:|---:|
| 100 | 0.523 | 0.057 | 0.521 | 0.054 |
| 500 | 0.989 | 0.045 | 0.988 | 0.051 |
| 1,000 | 1.000 | 0.053 | 1.000 | 0.051 |
| 5,000 | 1.000 | 0.051 | 1.000 | 0.052 |
| 10,000 | 1.000 | 0.043 | 1.000 | 0.050 |

Formula-specific empirical thresholds restored held-out null rejection to
approximately 0.05 for both methods, but the legacy threshold then grew with
$\sqrt n$. At the extreme `+I+G4` branch it rose from 3.11 at 100 sites to
18.54 at 10,000 sites; the refined threshold remained near 1.64.

After that independent calibration, the refined statistic also had modestly
higher power in the transition region. Representative paired results were:

| Model | $t$ | Sites | Legacy power | Refined power | Refined minus legacy, bootstrap 95% CI |
|---|---:|---:|---:|---:|---:|
| `+I+G4` | 50 | 1,000 | 0.315 | 0.365 | $0.050\ [0.041,0.059]$ |
| `+I+G4` | 50 | 5,000 | 0.823 | 0.885 | $0.062\ [0.054,0.070]$ |
| `+I+R4` | 20 | 500 | 0.391 | 0.438 | $0.047\ [0.039,0.056]$ |
| `+I+R4` | 20 | 1,000 | 0.646 | 0.682 | $0.036\ [0.028,0.044]$ |

The exact null-centered standardized separation was higher for the refined
score in 70 of 80 cells. The other ten cells were all at $t=10{,}000$, where
both finite-branch alternatives equal the null to numerical precision and
power is necessarily 0.05 after calibration. Thus the main improvement is
valid false-positive control; the power gain is real but smaller and localized
to the detection boundary.

#### Higher-precision figure rerun

The figure-sensitive conditions were rerun with 20,000 null and 20,000
alternative replicates per cell at

$$
t\in\{20,50,10000\}.
$$

This focused experiment contains 1.2 million paired simulated alignments. It
measures the shift directly as

$$
\Delta Z=Z_{\mathrm{refined}}-Z_{\mathrm{legacy}}
$$

and, after empirical calibration, as the difference in decision margins

$$
\Delta M
=
\left(Z_{\mathrm{refined}}-q_{0.95,\mathrm{refined}}\right)
-
\left(Z_{\mathrm{legacy}}-q_{0.95,\mathrm{legacy}}\right).
$$

At the extreme branch, the paired null shifts were:

| Sites | Mean $\Delta Z$, `+I+G4` | $\Delta$ false-positive rate | Mean $\Delta Z$, `+I+R4` | $\Delta$ false-positive rate |
|---:|---:|---:|---:|---:|
| 100 | -1.640 | -0.457 | -1.657 | -0.458 |
| 500 | -3.772 | -0.932 | -3.772 | -0.938 |
| 1,000 | -5.349 | -0.949 | -5.342 | -0.952 |
| 5,000 | -12.006 | -0.949 | -11.992 | -0.950 |
| 10,000 | -16.984 | -0.948 | -16.967 | -0.950 |

The negative $\Delta Z$ is desirable here: it removes the legacy positive-null
shift. It is not evidence of lost phylogenetic signal. Once each statistic is
measured relative to its independently trained null threshold, the refined
decision margin generally increases.

The focused calibrated-power results were:

| Model | $t$ | Sites | Legacy power | Refined power | Paired difference, bootstrap 95% CI |
|---|---:|---:|---:|---:|---:|
| `+I+G4` | 50 | 100 | 0.100 | 0.115 | $0.0148\ [0.0115,0.0182]$ |
| `+I+G4` | 50 | 1,000 | 0.317 | 0.348 | $0.0312\ [0.0261,0.0363]$ |
| `+I+G4` | 50 | 5,000 | 0.808 | 0.875 | $0.0673\ [0.0629,0.0715]$ |
| `+I+R4` | 20 | 100 | 0.151 | 0.164 | $0.0130\ [0.0093,0.0167]$ |
| `+I+R4` | 20 | 500 | 0.405 | 0.450 | $0.0442\ [0.0391,0.0488]$ |
| `+I+R4` | 20 | 1,000 | 0.634 | 0.684 | $0.0501\ [0.0454,0.0547]$ |

The apparent gain becomes small again when both methods approach either 5%
power or 100% power. The practical improvement is therefore concentrated
where the branch is difficult but still detectable.

#### Full dense branch-grid analysis

The manuscript run expands the focused analysis to 27 branch lengths:

$$
\begin{aligned}
t\in\{&
0.5,1,2,3,4,5,7.5,10,12.5,15,17.5,20,22.5,25,\\
&30,35,40,45,50,60,75,100,150,200,500,1000,10000
\}.
\end{aligned}
$$

The grid retains

$$
n\in\{100,500,1000,5000,10000\}
$$

and uses 20,000 null plus 20,000 alternative replicates in every cell. Across
the two models this is 270 model--branch--length cells and 10.8 million paired
simulated alignments. The 2.5-unit branch steps from 10 to 25 resolve the
+I+R4 transition, while the 5-unit steps from 30 to 50 and the additional
60/75 points resolve the +I+G4 transition.

The exact likelihood identity held with maximum error
$1.55\times10^{-11}$, and the refined exact null mean was below
$1.53\times10^{-16}$. Across all 270 cells:

| Quantity | Legacy counterfactual | Refined invariant null |
|---|---:|---:|
| Mean nominal null rejection | 0.9002 | 0.0492 |
| Nominal null-rejection range | 0.4998--1.0000 | 0.0423--0.0591 |
| Mean held-out calibrated rejection | 0.0501 | 0.0499 |

The dense steps reveal a diagonal power-gain ridge that the coarse grid
under-resolved. The largest calibrated gains were:

| Model | Sites | Branch | Legacy power | Refined power | Difference |
|---|---:|---:|---:|---:|---:|
| +I+G4 | 100 | 50 | 0.100 | 0.115 | 0.015 |
| +I+G4 | 500 | 40 | 0.342 | 0.384 | 0.042 |
| +I+G4 | 1,000 | 45 | 0.406 | 0.456 | 0.050 |
| +I+G4 | 5,000 | 60 | 0.482 | 0.584 | 0.102 |
| +I+G4 | 10,000 | 60 | 0.718 | 0.835 | 0.116 |
| +I+R4 | 100 | 10 | 0.558 | 0.573 | 0.015 |
| +I+R4 | 500 | 20 | 0.405 | 0.450 | 0.044 |
| +I+R4 | 1,000 | 22.5 | 0.449 | 0.522 | 0.073 |
| +I+R4 | 5,000 | 30 | 0.501 | 0.605 | 0.104 |
| +I+R4 | 10,000 | 35 | 0.402 | 0.565 | 0.163 |

The exact null-centered standardized separation was higher for the refined
score in 245 cells and equal in the 25 remaining extreme-saturation cells; it
was never lower. The largest negative empirical power difference was only
$-0.0059$ and occurred among 270 unadjusted comparisons, mainly near the null
floor. It was not supported by a lower exact separation.

### Why the relative-weighted curve drops abruptly

The relative weights are continuous functions of branch length, but the hard
MAP category partition is not. At branch length 30, the leading signed
contributions to pooled relative-weighted coherence were approximately

$$
+0.569\quad\text{from category 1},
\qquad
-0.187\quad\text{from category 2}.
$$

At branch length 50 they became

$$
+0.463\quad\text{and}\quad-0.399,
$$

which nearly cancel. Between branch lengths 75 and 100, category 2 becomes
empty and its sites move mainly into category 1; this is the largest adjacent
assignment change, with total-variation distance 0.407. Therefore the abrupt
pooled movement is caused by discrete category reassignment and cancellation,
not by an abrupt change in eigenvalue weights.

### Model ruler versus detection

For this G4 model, the model-based index was approximately 0.90 at branch
length 10, 0.98 at 30, and 0.996 at 50. Detectability nevertheless remained
high at 30 with 1,000 sites. This is not a contradiction: $S(t)$ measures the
fraction of model spectral energy lost, whereas the SatuTe tests ask whether
the remaining energy is detectable at the available sample size.

## Generated figures

The requested six-panel separation of effect, uncertainty, detection,
persistence, calibrated decisions and model saturation is:

![SatuTe saturation diagnostics](figures/generated/figure_saturation_diagnostics.png)

The hard-category mechanism is shown separately:

![Hard-category assignment and signed contributions](figures/generated/figure_saturation_category_assignment.png)

The alignment-length null verification is:

![Standard-normal null behavior across alignment lengths](figures/generated/figure_null_scaling_diagnostics.png)

The dedicated paired invariant-mixture grid is:

![Paired +I+G4 and +I+R4 correction grid](figures/generated/figure_weighted_invariant_grid.png)

The higher-precision direct-shift rerun is:

![Direct paired shift caused by the invariant correction](figures/generated/figure_weighted_invariant_shift.png)

The manuscript-scale dense branch-grid figure is:

![Full dense paired invariant-mixture analysis](figures/generated/figure_weighted_invariant_full.png)

The dense figure is exported as a vector
[PDF](figures/generated/figure_weighted_invariant_full.pdf) and
[SVG](figures/generated/figure_weighted_invariant_full.svg), plus a
4,800-by-4,200-pixel PNG.

The full manuscript run is reproduced with:

    python3 research/satute/experiments/003_invariant_mixture/run.py \
      --null-reps 20000 \
      --alternative-reps 20000 \
      --outdir research/satute/artifacts/local/weighted-invariant-full \
      --figure-prefix research/satute/artifacts/releases/generated-figures/figure_weighted_invariant_full_overview \
      --png-scale 3

    python3 research/satute/experiments/003_invariant_mixture/render_full.py --png-scale 3

## Reproduction

```bash
python3 research/satute/experiments/004_saturation_diagnostics/run.py \
  --iqtree ./build/iqtree3 \
  --sites 1000 \
  --reps 100 \
  --null-reps 500 \
  --workers 8
```

The invariant-mixture grid is reproduced with:

```bash
python3 research/satute/experiments/003_invariant_mixture/run.py
```

The higher-precision shift figure is reproduced with:

```bash
python3 research/satute/experiments/003_invariant_mixture/run.py \
  --branch-lengths 20,50,10000 \
  --null-reps 20000 \
  --alternative-reps 20000 \
  --outdir research/satute/artifacts/local/weighted-invariant-shift \
  --figure-prefix research/satute/artifacts/releases/generated-figures/figure_weighted_invariant_shift_overview

python3 research/satute/experiments/003_invariant_mixture/render_shift.py
```

The detailed and summarized data are written below
`research/satute/artifacts/local/saturation-diagnostics/`. The manuscript-ready
SVG, PDF and PNG are written as
`research/satute/artifacts/releases/generated-figures/figure_saturation_diagnostics.*`.
The category-mechanism figure is written as
`research/satute/artifacts/releases/generated-figures/figure_saturation_category_assignment.*`.
The null-scaling summary and figure are generated with
[`summarize_null_scaling.py`](summarize_null_scaling.py).

The exact verifier requires NumPy:

```bash
${PYTHON_BIN:-python3} research/satute/tests/exact/test_saturation_null.py
```
