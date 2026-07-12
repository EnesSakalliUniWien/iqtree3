# SatuTe Phase 1: Development Plan

Phase 1 includes:

- one fixed tree
- one alignment
- reversible stationary models only
- branch-level output
- testing all branches or a user-defined subset
- Bonferroni correction
- discrete rate-category models such as `+G4`
- three formula rows: `dominant`, `eigenvalue_weighted`, and
  `mixture_likelihood_weighted`

For the current code-facing implementation status, assumptions, verification
coverage and limitations, see
[port-status.md](port-status.md).
For a requirement-by-requirement audit of what is proven versus still missing,
see
[completion-audit.md](completion-audit.md).

For each tested branch, compute a SatuTe summary that answers:

> Is there still detectable phylogenetic signal linking the two subtrees separated by this branch?

The native implementation produces:

- estimated coherence
- estimated standard error
- z-score
- `alpha`
- `alpha_adjust`
- classification: `informative` or `saturated`
- the formula name
- the IQ-TREE rate category, IQ-TREE rate multiplier, and effective branch length
- spectral modes, eigenvalues, and formula weights
- the formula-independent model information fraction and monotonic saturation
  index

### 1. Branch-side partial likelihoods

The current phase-1 implementation computes ordinary scalar partial likelihood
vectors directly for the two sides of every tested branch. This deliberately
avoids reading IQ-TREE's packed/vectorized `partial_lh` buffers as if they were
site-contiguous scalar arrays.

For each branch side and site pattern, the implementation:

- prunes the subtree away from the focal branch in ordinary state space
- builds transition matrices from the fitted model eigensystem
- converts the side likelihood vector into posterior state probabilities with
  the model stationary frequencies

### 2. SatuTe Formula

For a tested branch and a site pattern, let `L_left(s)` and `L_right(s)` be the
ordinary partial likelihood vectors for the two branch sides, and let `pi(s)` be
the stationary state frequency. The side posteriors are:

```text
P_left(s)  = L_left(s)  * pi(s) / sum_x L_left(x)  * pi(x)
P_right(s) = L_right(s) * pi(s) / sum_x L_right(x) * pi(x)
```

For every dominant non-zero eigenmode `k`, the side factors are:

```text
F_left(k)  = sum_s right_eigenvector(s,k) * P_left(s)
F_right(k) = sum_s right_eigenvector(s,k) * P_right(s)
```

The site coherence is:

```text
C_site = sum_k F_left(k) * F_right(k)
```

The reported `satC` is the frequency-weighted mean over valid site patterns.
The null variance uses the product of the left and right second moments over
the same modes. The `dominant` formula uses the published SatuTe dominant
eigenspace. The `eigenvalue_weighted` formula uses every non-stationary mode with relative
weight `exp((lambda_k - lambda_1) * t)`, where `t` is the tested branch length
after any rate-category rescaling.

The independent comparison script also computes an offline `all_unweighted`
formula with every non-stationary mode assigned weight one. It is not part of
the native `.sat.stat` interface.

The refined `mixture_likelihood_weighted` formula avoids hard rate-category
assignment. For positive-rate category prior `q_c`, category-specific
independent-side likelihoods `A_sc` and `B_sc`, and site/category signal
`C_sc`, it uses

$$
\tau_{sc}=\frac{q_c A_{sc}B_{sc}}{\sum_d q_d A_{sd}B_{sd}},
\qquad
C_s=\sum_c\tau_{sc}C_{sc},
$$

with spectral persistence weights

$$
w_{kc}=\exp(\lambda_k r_c t).
$$

The native calculation subtracts one common log-weight shift from every
$\lambda_kr_ct$ before exponentiation. This prevents underflow while preserving
the relative weights across categories and leaves the standardized statistic
unchanged. For an explicit zero-rate invariant category, the null denominator
uses the invariant joint likelihood

$$
J_{s0}=\sum_a\pi_aL_{A,0}(a)L_{B,0}(a),
$$

and the invariant component contributes no numerator signal because its focal
transition is identical at finite and infinite branch lengths. The variance is
the ordinary frequency-weighted sample variance of the final site values
$C_s$, so category uncertainty and between-category variation are retained.

### 3. Rate categories

For models with discrete rate heterogeneity, native `--satute` does not estimate
rate categories itself. It uses IQ-TREE's own empirical-Bayes assignment through
`computePatternRates`: each site pattern is assigned to its maximum-posterior
category, and each category uses the fitted IQ-TREE multiplier returned by
`getRate(c)`. This is the same category/rate information exposed through
IQ-TREE's `.rate` output. Pure `+I` is deliberately not converted into
category-wise SatuTe output because IQ-TREE's own `.rate` writer emits no
per-site category rows for a pure invariable-sites model. Models such as `+G`,
`+R`, and IQ-TREE-supported `+I+G`/`+I+R` combinations are handled when
`computePatternRates` provides pattern categories. The original Python SatuTe
workflow obtained the same information indirectly by running IQ-TREE with
`-wspr`, reading the `.siteprob` posterior columns and assigning sites to the
largest posterior category.

For `dominant` and `eigenvalue_weighted`, each rate category is evaluated on a
tree whose branch lengths are multiplied by the IQ-TREE category rate. The
`.sat.stat` file contains one pooled row and one row per category for each branch
and formula. Their pooled row uses the SatuTe global-variance rule:

```text
satC_pooled   = sum_c n_c / n * satC_c
satVar_pooled = sum_c n_c / n * satVar_c
satSE_pooled  = sqrt(satVar_pooled / n)
```

`mixture_likelihood_weighted` instead emits one inherently pooled row, using
all IQ-TREE category rates and proportions at every site through the soft
responsibilities above.

Every row also reports the model-based scale

$$
\operatorname{InformationFraction}(t)
=\sum_c q_c\frac{1}{K}\sum_{k\ne0}\exp(2\lambda_k r_ct),
\qquad
\operatorname{SaturationIndex}(t)=1-\operatorname{InformationFraction}(t),
$$

with category proportions normalized if necessary. Category rows use their own
effective length; pooled rows integrate over the IQ-TREE proportions. These
columns are identical across formulas because they describe the fitted model,
not observed test power.

For homogeneous models, only pooled rows are emitted.

### 4. Reversible model eigensystem

Use:

- eigensystem decomposition in [modelmarkov.h](../../model/modelmarkov.h)
- eigenvalue/eigenvector accessors in [modelmarkov.h](../../model/modelmarkov.h)

### 5. Reversible-kernel gate

Use:

- [modelsubst.h](../../model/modelsubst.h)

Phase 1 rejects unsupported cases early: supertrees, mixture models and
ascertainment-bias correction.

### 6. Simulated examples

Run the full phase-1 verification bundle with:

```bash
research/satute/tests/native/run_satute_phase1_checks.sh ./build/iqtree3
```

This runs the native smoke test, the IQ-TREE-owned rate-category regression, the
formula-comparison reference check, and a larger-tree reference check into one
output root. It also exercises the standalone Python reference CLI on one
homogeneous GTR case, one GTR+G4 `.rate` category case, and one six-taxon GTR+G4
case, then compares those outputs back to native `.sat.stat`.
The formula comparison fails if native `dominant` or `eigenvalue_weighted`
pooled rows diverge from the independent Python reference
for `satC`, `satVar`, `satSE`, `satZ`, `satP`, the decision or the number of
modes.

Run the reusable Python reference directly with:

```bash
python3 research/satute/tests/native/satute_reference.py \
  --alignment /path/to/alignment.fa \
  --tree /path/to/tree-or-sat.tree \
  --model 'GTR{1,2,1,1,2,1}+F{0.30,0.20,0.20,0.30}' \
  --split A,B \
  --compare-sat-stat /path/to/native.sat.stat
```

When native `.sat.stat` output is available, an exact IQ-TREE branch can also be
selected by ID:

```bash
python3 research/satute/tests/native/satute_reference.py \
  --alignment /path/to/alignment.fa \
  --tree /path/to/tree-or-sat.tree \
  --model 'GTR{1,2,1,1,2,1}+F{0.30,0.20,0.20,0.30}+G4{0.5}' \
  --rate-file /path/to/native.rate \
  --iqtree-report /path/to/native.iqtree \
  --branch-id 7 \
  --compare-sat-stat /path/to/native.sat.stat
```

The reference CLI now parses a general fixed Newick tree, finds the requested
branch by split, and recursively recomputes ordinary side partials away from
that edge. It supports nucleotide homogeneous JC/GTR rows and optional IQ-TREE
`.rate` category rows via `--rate-file`; the soft-mixture reference also reads
category rates and proportions from `--iqtree-report`. With
`--compare-sat-stat`, the same CLI
fails if the recomputed Python reference rows disagree with native `.sat.stat`.
The `--split` argument may name either side of the same branch, for example
`A,B` or `C,D`; non-matching splits are rejected. The `--branch-id` argument
uses native `.sat.stat` to recover the exact split for one IQ-TREE branch ID, and
then compares back to rows for that same ID. The phase-1 bundle checks this by
requiring `--split A,C` to fail on an `A,B|C,D` tree and by matching a six-taxon
GTR+G4 split and branch ID against native output.

Run the JC and GTR smoke examples with:

```bash
research/satute/tests/native/test_satute_simulated.sh ./build/iqtree3
```

This smoke test also reruns one branch through `--satute-edges` and
`--satute-alpha` to verify branch-subset output and the requested significance
level in `.sat.stat`. It also requires an absent branch ID to fail with
`Requested SatuTe branch IDs not found: 999999`.

Run the IQ-TREE-owned rate-category regression with:

```bash
research/satute/tests/native/test_satute_rate_categories.py ./build/iqtree3
```

This check simulates one alignment, runs native `--satute --rate` under `+G`,
`+R`, `+I+G`, and `+I+R`, and verifies that nonzero category counts and
category rates in `.sat.stat` match IQ-TREE's `.rate` output. For `JC+G4`, it
also recomputes the `A,B` branch dominant SatuTe category rows and pooled row
with an independent Python reference calculation. For `GTR+G4`, it recomputes
the same branch for `dominant` and `eigenvalue_weighted`. It also checks
that pure `+I` fails with an explicit message because IQ-TREE emits no per-site
category rows for that model.

The GTR examples use:

```text
GTR{1,2,1,1,2,1}+F{0.30,0.20,0.20,0.30}
```

For this setup, a focal branch length of `0.50` should be called informative,
whereas `8.00` should be called saturated.

For a comparison against an all-nonzero-mode eigenvector variant, see
[formula-comparison.md](formula-comparison.md).

Native smoke checks after the rate-category port used:

```bash
build/iqtree3 -s /tmp/satute-native-smoke/jc.fa \
  -te /tmp/satute-native-smoke/tree.nwk \
  -m JC --prefix /tmp/satute-native-smoke/jc_sat \
  -T 1 --redo --quiet --satute

build/iqtree3 -s /tmp/satute-native-smoke/gtrg4.fa \
  -te /tmp/satute-native-smoke/tree.nwk \
  -m 'GTR{0.6676,3.7807,4.2833,0.5354,0.8718,1.0}+F{0.125,0.436,0.191,0.245}+G4{0.5}' \
  --prefix /tmp/satute-native-smoke/gtrg4_sat \
  -T 1 --redo --quiet --satute
```

The `JC` run emits identical `dominant` and `eigenvalue_weighted`
pooled rows. The `GTR+G4` run emits pooled rows plus four category rows per
formula, and the eigenvalue weights differ by IQ-TREE category rate.

## Key Implementation Fact

For reversible kernels, IQ-TREE stores internal branch partials in vectorized
eigen-space buffers.

Relevant code:

- inverse-eigenvector projection in [phylokernelnew.h](../../tree/phylokernelnew.h)
- branch likelihood assembly in [phylokernelnew.h](../../tree/phylokernelnew.h)

Those buffers are not scalar site-contiguous, so phase 1 does not read them
directly for the SatuTe statistic. It recomputes the ordinary side partials in
scalar state space and then applies the posterior/eigenvector formula above.

## What Data Is Needed For One Branch

For a tested branch `AB`, use:

- the two directed edge objects for that branch
- ordinary partial likelihood vectors recomputed from side A and side B
- the model eigenvalues
- the model right eigenvectors
- the model stationary frequencies
- pattern frequencies
- the number of taxa on each side of the branch

Where they come from:

- edge partials: [phylonode.h](../../tree/phylonode.h)
- branch likelihood trigger: [phylotree.h](../../tree/phylotree.h)
- eigenvalues/eigenvectors: [modelmarkov.h](../../model/modelmarkov.h)
- pattern/category buffers: [phylotree.h](../../tree/phylotree.h)

Do not use these as the core statistic:

- bootstrap values
- concordance factors
- ancestral-state reconstruction
- overall tree log-likelihood alone
- a new custom traversal
