# SatuTe And The Eigenvalue-Weighted Extension

Saturation in phylogenetics describes the loss of historical signal after
repeated substitutions. SatuTe makes this question branch-specific. For a
branch that separates two subtrees, the method tests whether the site patterns
on the two sides still show dependence beyond what is expected under
independent evolution. This formulation allows one branch to be saturated while
other parts of the same alignment remain informative.

For one site, each side of the focal branch produces a partial likelihood vector
over possible endpoint states. In a DNA model, this vector has four entries.
SatuTe converts each branch-side likelihood vector into a posterior state
distribution by multiplying by the stationary frequencies and normalizing. The
two sides can then be compared in the same model-defined state space.

The substitution model supplies the spectral coordinates for this comparison.
The stationary eigenmode represents equilibrium and is not used as a
phylogenetic contrast. Non-stationary eigenvectors define directions of state
contrast, and their eigenvalues define the expected rate of decay along a
branch. SatuTe projects the left and right posterior state distributions onto
the relevant non-stationary eigenvectors and multiplies the paired projections.
Under subtree independence, the expected value of this product is zero.

In the JC model, all three non-stationary eigenvalues are equal. The dominant
non-stationary eigenspace therefore contains all three directions, and the
published statistic and the eigenvalue-weighted statistic give the same test
after normalization. This provides a useful reference case: if the model has no
spectral separation among non-stationary directions, eigenvalue weighting
should not create one.

In a typical GTR model, the non-stationary eigenvalues differ. The published
SatuTe statistic uses the slowest-decaying non-stationary direction. Other
eigenvectors may still carry agreement between the two subtrees, but their
eigenvalues indicate faster decay. The eigenvalue-weighted statistic includes
those directions with weights determined by their expected persistence across
the focal branch.

The proposed site-level weighted statistic has the form
`sum_k exp((lambda_k - lambda_*) * branch_length) * C_k`, where `lambda_*` is
the dominant non-stationary eigenvalue and `C_k` is the coherence coefficient
for mode `k`. The relative weight is one for the dominant mode and less than
one for faster-decaying modes. This keeps the published SatuTe coefficient as
the reference scale while allowing additional modes to contribute when the
branch is short enough for those modes to remain relevant under the fitted
model.

The simulation analysis must preserve the structure of the original SatuTe
study. The Fig. 2 design used JC simulations on a five-taxon external branch
and a 16-taxon internal branch, varied the focal branch over sixteen values,
used site lengths of 100, 1000 and 10000, and generated 1000 alignments per
parameter combination. Each alignment was analyzed on the true tree with fixed
lengths, on the true topology with maximum-likelihood lengths, on a
maximum-likelihood inferred tree and on the same inferred tree with the
Bonferroni correction.

The reduced runs in this repository test the implementation, target-branch
selection and output parsing. They are not estimates of power or calibration.
The method comparison should be made on paired replicates: the same simulated
alignment, target branch, inferred or supplied tree and significance threshold
should be used for the dominant and eigenvalue-weighted statistics.

At the current stage, the published dominant statistic remains the reference
method. The eigenvalue-weighted statistic is the extension to evaluate because
it uses both eigenvectors and eigenvalues. Its status depends on the full
simulation analysis, including type I error, power and behavior under
maximum-likelihood tree selection.

An explicit invariant-site component requires an additional correction. A
zero-rate category retains the same ancestral state across the focal branch
even in the infinite-branch limit. Its null contribution is therefore a joint
subtree likelihood rather than the product of the two subtree likelihoods. The
refined soft-mixture statistic places this invariant joint term in the null
denominator and allows only positive-rate categories to contribute
finite-versus-infinite-branch signal.

The full paired +I analysis used +I+G4 and +I+R4 models, five alignment
lengths, 27 branch lengths and 20,000 null plus 20,000 alternative replicates
per cell. The refined nominal false-positive rate averaged 4.92%, whereas the
counterfactual invariant construction averaged 90.0%. After independent null
calibration, the largest refined power gains were 11.6 percentage points for
+I+G4 and 16.3 percentage points for +I+R4. These gains form a narrow ridge at
the detection boundary; the primary benefit is restoration of the correct
null.
