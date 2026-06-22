# Formula From Scratch

Consider a branch that splits a tree into two subtrees. Call the two sides
`left` and `right`. For one site pattern, pruning each side away from the focal
branch gives two ordinary partial likelihood vectors:

```text
L_left(s)
L_right(s)
```

where `s` is a nucleotide state. These vectors say how compatible each possible
state at the branch endpoint is with the observed descendants on that side.

SatuTe converts these likelihood vectors into posterior state distributions by
using the stationary distribution of the substitution model:

```text
P_left(s)  = L_left(s)  * pi(s) / sum_x L_left(x)  * pi(x)
P_right(s) = L_right(s) * pi(s) / sum_x L_right(x) * pi(x)
```

The posterior step matters because the branch-side likelihood vector is not yet
a state distribution. Multiplication by `pi` makes the vector interpretable as
the posterior probability of the hidden branch-end state under the model.

The substitution model decomposes state-space variation into eigenmodes. The
stationary eigenmode is the equilibrium component and does not carry
phylogenetic contrast. The non-stationary eigenvectors define directions of
state contrast, and the associated eigenvalues define how quickly each
direction decays along a branch.

For each non-stationary mode `k`, define:

```text
F_left(k)  = sum_s right_eigenvector(s,k) * P_left(s)
F_right(k) = sum_s right_eigenvector(s,k) * P_right(s)
```

The published SatuTe statistic uses the dominant non-stationary eigenspace,
meaning the non-stationary eigenvalue closest to zero. This is the
slowest-decaying signal direction. For one site:

```text
C_site = sum_{k in dominant eigenspace} F_left(k) * F_right(k)
```

The reported coherence `satC` is the average of `C_site` over site patterns,
weighted by pattern frequency. Under subtree independence, the expected
coherence is zero. The variance estimator is constructed from the second
moments of the left and right factors, which gives the standard error and the
one-sided z-test.

## JC Versus GTR

Under JC, the three non-stationary eigenvalues are equal. The dominant
non-stationary eigenspace has multiplicity three, so all three non-stationary
eigenvectors are used. Any all-mode variant collapses back to the same test
because there are no faster or slower non-stationary modes to separate.

Under a typical GTR model, the eigenvalues split. Usually only one
non-stationary eigenvalue is closest to zero. The published SatuTe statistic
therefore uses one dominant direction. Other eigenvectors still describe real
state contrasts, but their more negative eigenvalues mean that those contrasts
decay faster along branches.

This distinction motivates the eigenvalue-weighted diagnostic. If non-dominant
GTR eigenvectors are included, their decay rates need to enter the statistic.
