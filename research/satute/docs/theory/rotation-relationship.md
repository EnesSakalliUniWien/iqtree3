# Relationship to *When Reversible Mutation Becomes Rotation*

## Result

The rotation project and SatuTe use the same reversible spectral backbone, but
they use it for different tasks:

- the rotation project is a forward construction from a rate matrix to exact
  parent-conditioned child probabilities and a branch circuit;
- SatuTe is an inverse diagnostic of whether an observed alignment retains
  detectable coherence across a fitted branch;
- the new `SaturationIndex` is neither of those tests. It is a deterministic,
  monotonic coordinate of the fitted branch model.

This separation is the main lesson from comparing the projects. The model's
loss of spectral persistence, the alignment's remaining coherence and the
standardized evidence for that coherence should not be represented by one
number.

## Common spectral object

Let a reversible, irreducible rate matrix have stationary distribution
$\boldsymbol\pi$ and symmetric representation

$$
\widehat Q
=
\Pi^{1/2}Q\Pi^{-1/2}
=
\widehat S
\operatorname{diag}(0,\lambda_2,\ldots,\lambda_{K+1})
\widehat S^\top,
$$

where $\lambda_k<0$ are the $K$ nonstationary eigenvalues. For rate category
$c$, define the modal persistence

$$
\rho_{kc}(t)=\exp(\lambda_k r_ct).
$$

The symmetrized transition operator is

$$
\widehat P_c(t)
=
\widehat S
\operatorname{diag}(1,\rho_{2c}(t),\ldots,\rho_{K+1,c}(t))
\widehat S^\top.
$$

Writing

$$
\Pi_0=\sqrt{\boldsymbol\pi}\sqrt{\boldsymbol\pi}^{\,\top}
$$

for its stationary component gives the basis-independent identity

$$
\frac{1}{K}
\left\|\widehat P_c(t)-\Pi_0\right\|_F^2
=
\frac{1}{K}\sum_{k\ne0}\rho_{kc}(t)^2.
$$

Therefore the native information fraction and saturation index are

$$
R(t)
=
1-S(t)
=
\frac{1}{\sum_cq_c}
\sum_cq_c\frac{1}{K}\sum_{k\ne0}
\exp(2\lambda_kr_ct),
$$

$$
S(t)=1-R(t).
$$

This proves that the choice $a_k=1$ is not arbitrary: it is the normalized
squared Frobenius energy of the nonstationary part of the reversible transition
operator. It is invariant to eigenvector sign choices and to rotations within
a repeated-eigenvalue subspace.

For a rate mixture, this implementation averages the category-specific squared
operator distances,

$$
\operatorname E_c
\left[\left\|\widehat P_c(t)-\Pi_0\right\|_F^2\right],
$$

not the squared distance of the category-averaged transition operator. Those
are different quantities. The implemented version answers how much spectral
energy remains for a site drawn from the fitted rate distribution.

## Rotation interpretation

The manuscript defines spectral phase angles by

$$
\rho_{kc}(t)=\exp(\lambda_kr_ct)=\cos\theta_{kc}(t),
\qquad 0\leq\theta_{kc}(t)\leq\frac{\pi}{2}.
$$

Consequently,

$$
R(t)
=
\frac{1}{\sum_cq_c}
\sum_cq_c\frac{1}{K}\sum_{k\ne0}
\cos^2\theta_{kc}(t),
$$

and

$$
\boxed{
S(t)
=
\frac{1}{\sum_cq_c}
\sum_cq_c\frac{1}{K}\sum_{k\ne0}
\sin^2\theta_{kc}(t)
}.
$$

Thus $S(t)$ is the mean squared spectral rotation away from the initial
nonstationary axes, while $R(t)$ is their mean squared persistence.

For JC2, where the only nonstationary eigenvalue is $-2\mu$,

$$
S(t)=1-e^{-4\mu t}=\sin^2\theta(t).
$$

For homogeneous JC69 in IQ-TREE's expected-substitution normalization,
$\lambda_2=\lambda_3=\lambda_4=-4/3$, so

$$
S(t)=1-e^{-8t/3}.
$$

The native closed-form regression test checks this identity.

### Important boundary

The spectral angles $\theta_k(t)$ are not the same as the compiled
parent-conditioned entangling angle. In the generic four-state construction,
measuring the spectral unitary would produce squared matrix-entry magnitudes,
not $P_{ab}(t)$. The rotation manuscript therefore continues through

$$
Q
\longrightarrow
F_a\delta(t)
\longrightarrow
v_a(t)
\longrightarrow
C_Q(t),
$$

where $\delta_k(t)=1-\rho_k(t)$ and
$|v_a(t)_b|^2=P_{ab}(t)$. The interpretation of $S(t)$ is an exact statement
about the spectral phase coordinates and transition-operator energy. It is not
a claim that $S(t)$ is a child-state probability, an entanglement measure or a
quantum-advantage metric.

## Map of the current methods

| Method | Spectral quantity | What it answers | Why it need not move monotonically |
|---|---|---|---|
| Published `dominant` | Coherence in the slowest-decaying eigenspace | Is the most persistent observed contrast detectable across the branch? | The observed signed coherence and its estimated standard error both vary with the alignment. |
| `eigenvalue_weighted` | $w_{kc}=\rho_{kc}/\rho_{*c}=e^{(\lambda_k-\lambda_*)r_ct}$ | Do all modes, weighted by persistence relative to the slowest mode, provide detectable coherence? | Dividing by $\rho_{*c}$ removes the common radial decay. Hard rate-category reassignment and cancellation of signed category contributions can dominate the pooled curve. |
| `SaturationIndex` | $1-\operatorname E_c[K^{-1}\sum_k\rho_{kc}^2]$ | How far has the fitted branch model progressed from zero-length persistence toward its stationary limit? | It is monotonic under the fixed-model assumptions; it is not calculated from observed coherence. |

The relative-weighted vector

$$
\left(\frac{\rho_{2c}}{\rho_{*c}},\ldots,
      \frac{\rho_{K+1,c}}{\rho_{*c}}\right)
$$

describes the *shape* of the surviving spectrum after its common scale has been
removed. The new information fraction $R(t)$ supplies the missing radial
magnitude. This explains why the relative-weighted statistic can still look
stable or move abruptly after the branch is already close to model saturation:
it was deliberately designed as a normalized detection statistic, not as a
radial saturation coordinate.

## What should be retained from each project, in priority order

1. **Keep model saturation separate from statistical detection.** Report
   `SaturationIndex` beside `satC`, `satSE`, `satZ` and calibrated decisions.
   Do not turn `satZ` into a saturation ruler.

2. **Use $S(t)$ as the common branch coordinate across models.** A simulation
   or benchmark grid matched on $S(t)$ compares the same fraction of remaining
   spectral energy more fairly than a grid matched only on raw branch length.
   Rotation-circuit accuracy, SatuTe power and ancestral reconstruction can all
   be plotted against this common coordinate.

3. **Retain the full modal spectrum when anisotropy matters.** The scalar
   $S(t)$ is a radius. It cannot distinguish spectra with the same squared
   energy but different slow/fast-mode composition. Report the persistence
   vector $\boldsymbol\rho(t)$ or an explicit anisotropy summary beside $S(t)$
   when comparing GTR models.

4. **Treat rate mixing at the level of the scientific target.** Averaging
   squared category persistence is appropriate for expected information in a
   randomly drawn site category. Squaring a category-averaged transition would
   instead measure the distance of the marginal transition law and includes
   cross-category terms. The two definitions should not be interchanged.

5. **Preserve eigenspace invariance.** Equal coefficients within a repeated
   eigenvalue block make $S(t)$ independent of the numerical eigenbasis. Any
   future mode-specific coefficients $a_k$ need to be constant within
   degenerate eigenspaces or defined through spectral projectors.

6. **Use the rotation construction as a forward verification oracle, not as a
   replacement test.** Its exact $Q\to P(t)\to C_Q(t)$ checks can validate the
   same transition law used by SatuTe. The alignment-based null, variance and
   power questions still require the SatuTe simulation suite.

## Recommended analysis layout

For every branch or simulation condition, plot:

1. $S(t)$ or $R(t)=1-S(t)$;
2. all category-specific $\rho_{*c}(t)$ values;
3. raw coherence $\bar C(t)$ for each detection method;
4. $\operatorname{SE}\{\bar C(t)\}$;
5. $Z(t)$;
6. calibrated probability of an informative decision;
7. the modal persistence vector or an anisotropy summary for non-JC models.

The first, second and seventh panels describe the fitted model. The remaining
panels describe the observed statistic and its sampling behavior.

## Sources reviewed in the rotation repository

The relationship above was checked against the current manuscript and its
implementation-facing documentation:

- [manuscript abstract and structure](<../../../Is it rotation/manuscript/main.tex>)
- [reversible-rate introduction](<../../../Is it rotation/manuscript/sections/01_introduction/index.tex>)
- [JC2 decay and rotation](<../../../Is it rotation/manuscript/sections/03_jc2_two_state_model/index.tex>)
- [GTR symmetrization](<../../../Is it rotation/manuscript/sections/05_four_state_models/gtr/01_model_and_symmetrization.tex>)
- [four-state branch construction](<../../../Is it rotation/manuscript/sections/06_four_state_gate_construction/index.tex>)
- [spectral-unitary versus sampler distinction](<../../../Is it rotation/manuscript/sections/06_four_state_gate_construction/01_spectral_gate_map.tex>)
- [results](<../../../Is it rotation/manuscript/sections/09_results/index.tex>)
- [discussion](<../../../Is it rotation/manuscript/sections/10_conclusion/index.tex>)
- [compiled manuscript](<../../../Is it rotation/build/main.pdf>)

No files in the rotation repository were changed during this review.
