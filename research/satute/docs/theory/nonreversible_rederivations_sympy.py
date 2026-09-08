#!/usr/bin/env python3
"""Symbolic verification of the three non-reversible re-derivations.

Companion to rotation-relationship.md, section
"Non-reversible extension: three re-derivations". Every claim of that
section is rebuilt here in sympy and asserted, exactly where the algebra
is exact and to 30 digits where radicals appear.

Run with an environment that provides sympy, for example the rotation
repository environment:

    "../../../../../Is it rotation/.venv/bin/python" \
        nonreversible_rederivations_sympy.py

The worked example is a fixed irreducible, non-reversible three-state
rate matrix with rational entries, chosen so that the two nonstationary
eigenvalues form the complex conjugate pair -6 +/- i*sqrt(2): a genuine
damped rotation.
"""

import sympy as sp

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    print(f"PASS  {name}")


def is_zero(expr, digits=30, tol=sp.Float("1e-25")):
    """Exact zero if sympy can prove it, else a 30-digit numerical zero."""
    simplified = sp.simplify(expr)
    if simplified == 0:
        return True
    return abs(sp.Abs(simplified).evalf(digits)) < tol


def zero_matrix(matrix):
    return all(is_zero(entry) for entry in matrix)


# --------------------------------------------------------------------------
# The worked example: irreducible, non-reversible, rational
# --------------------------------------------------------------------------

n = 3
Q = sp.Matrix([
    [-5, 4, 1],
    [1, -3, 2],
    [3, 1, -4],
])

check("rows of Q sum to zero", Q * sp.ones(n, 1) == sp.zeros(n, 1))

# Stationary distribution: the normalized left null vector of Q.
pi_vec = Q.T.nullspace()[0]
pi_vec = pi_vec / sum(pi_vec)
check("pi is rational and positive", all(p > 0 for p in pi_vec))
check("pi^T Q = 0", (pi_vec.T * Q) == sp.zeros(1, n))

Pi = sp.diag(*pi_vec)
Pi12 = sp.diag(*[sp.sqrt(p) for p in pi_vec])
Pim12 = Pi12.inv()
sqrt_pi = sp.Matrix([sp.sqrt(p) for p in pi_vec])
Pi0 = sqrt_pi * sqrt_pi.T

# Symmetrized generator and its entrywise definition.
Qh = Pi12 * Q * Pim12
Qh_entrywise = sp.Matrix(
    n, n, lambda x, y: sp.sqrt(pi_vec[x] / pi_vec[y]) * Q[x, y]
)
check("Qh entries are sqrt(pi_x/pi_y) Q_xy", zero_matrix(Qh - Qh_entrywise))

check("sqrt(pi) is a right null vector of Qh", zero_matrix(Qh * sqrt_pi))
check("sqrt(pi) is a left null vector of Qh", zero_matrix(sqrt_pi.T * Qh))

# --------------------------------------------------------------------------
# Derivation 1: antisymmetric part and damped rotations
# --------------------------------------------------------------------------
print("\n-- 1. Antisymmetric part and damped rotations --")

A_hat = (Qh + Qh.T) / 2
W_hat = (Qh - Qh.T) / 2

# Time reversal and the additive reversibilization.
Q_star = Pi.inv() * Q.T * Pi
check(
    "A_hat = Pi^(1/2) (Q + Q*)/2 Pi^(-1/2)",
    zero_matrix(A_hat - Pi12 * (Q + Q_star) / 2 * Pim12),
)
check(
    "W_hat = Pi^(1/2) (Q - Q*)/2 Pi^(-1/2)",
    zero_matrix(W_hat - Pi12 * (Q - Q_star) / 2 * Pim12),
)

# Detailed balance pi_x Q_xy = pi_y Q_yx is exactly W_hat = 0.
detailed_balance = sp.Matrix(
    n, n, lambda x, y: pi_vec[x] * Q[x, y] - pi_vec[y] * Q[y, x]
)
check("this Q is non-reversible", not zero_matrix(detailed_balance))
check("W_hat is nonzero for this Q", not zero_matrix(W_hat))
reversibilization = (Q + Q_star) / 2
db_rev = sp.Matrix(
    n, n,
    lambda x, y: pi_vec[x] * reversibilization[x, y]
    - pi_vec[y] * reversibilization[y, x],
)
check("(Q + Q*)/2 satisfies detailed balance", zero_matrix(db_rev))

# Complex conjugate eigenvalue pair with negative real part.
eigenvalues = list(Q.eigenvals())
nonzero_eigenvalues = [lam for lam in eigenvalues if sp.simplify(lam) != 0]
check("two nonstationary eigenvalues", len(nonzero_eigenvalues) == 2)
check(
    "eigenvalues are the conjugate pair -6 +/- i sqrt(2)",
    {sp.simplify(lam) for lam in nonzero_eigenvalues}
    == {-6 + sp.I * sp.sqrt(2), -6 - sp.I * sp.sqrt(2)},
)

# tr(N^T W_hat N) = 0 for a fully generic N: the rotation never dissipates.
N_generic = sp.Matrix(n, n, sp.symbols("n0:9"))
check(
    "tr(N^T W_hat N) = 0 for generic N",
    sp.expand((N_generic.T * W_hat * N_generic).trace()) == 0,
)
check(
    "tr(N^T Qh N) = tr(N^T A_hat N) for generic N",
    sp.expand(
        (N_generic.T * Qh * N_generic).trace()
        - (N_generic.T * A_hat * N_generic).trace()
    )
    == 0,
)
# Together with dN/dt = Qh N (which holds because Qh Pi0 = 0), this gives
# d/dt ||Ph(t) - Pi0||_F^2 = 2 tr(N^T A_hat N) <= 0.
check("Qh Pi0 = 0", zero_matrix(Qh * Pi0))


def expm_numeric(matrix, terms=80):
    """Truncated exponential series; ample precision for the norms here."""
    result = sp.eye(n).evalf(30)
    term = sp.eye(n).evalf(30)
    for k in range(1, terms + 1):
        term = term * matrix / k
        result = result + term
    return result


def frobenius_energy(t_value):
    difference = expm_numeric(Qh.evalf(30) * t_value) - Pi0.evalf(30)
    return sum(entry**2 for entry in difference)


times = [sp.Rational(1, 5), sp.Rational(2, 5), sp.Rational(3, 5), sp.Rational(4, 5)]
energies = [frobenius_energy(t) for t in times]
check(
    "Frobenius energy decreases monotonically",
    all(energies[i] > energies[i + 1] for i in range(len(energies) - 1)),
)

# The eigenvalue-sum formula fails: sum_k |rho_k(t)|^2 = 2 exp(-12 t) here,
# and Qh is not normal, so the Frobenius energy differs from it.
t_probe = sp.Rational(3, 10)
eigen_sum = (2 * sp.exp(-12 * t_probe)).evalf(30)
energy_probe = frobenius_energy(t_probe)
check(
    "Qh is not normal",
    not zero_matrix(Qh * Qh.T - Qh.T * Qh),
)
check(
    "Frobenius energy differs from the eigenvalue sum",
    abs(energy_probe - eigen_sum) / eigen_sum > sp.Float("1e-9"),
)

# --------------------------------------------------------------------------
# Derivation 2: directed edges and canonical correlations
# --------------------------------------------------------------------------
print("\n-- 2. Directed edges and canonical correlations --")

# The pulley principle is the identity Q^T Pi = Pi Q, equivalently Q = Q*.
check("Q^T Pi = Pi Q* holds always", zero_matrix(Q.T * Pi - Pi * Q_star))
check(
    "Q^T Pi != Pi Q, so J(t1,t2) leaves Pi P(t1+t2) at first order",
    not zero_matrix(Q.T * Pi - Pi * Q),
)

# M(t1,t2) = Pi^(-1/2) J Pi^(-1/2) = Ph(t1)^T Ph(t2): a generic identity.
generic_P1 = sp.Matrix(n, n, sp.symbols("p0:9"))
generic_P2 = sp.Matrix(n, n, sp.symbols("q0:9"))
s = sp.symbols("s1:4", positive=True)
D12 = sp.diag(*s)          # a generic Pi^(1/2)
D = D12**2                 # the corresponding Pi
Dm12 = D12.inv()
check(
    "M = Ph(t1)^T Ph(t2) for generic P(t1), P(t2), Pi",
    sp.expand(
        Dm12 * generic_P1.T * D * generic_P2 * Dm12
        - (D12 * generic_P1 * Dm12).T * (D12 * generic_P2 * Dm12)
    )
    == sp.zeros(n, n),
)

# sqrt(pi) is a left and right singular vector of M with sigma_0 = 1.
# The order-2 series of P(t) preserves P 1 = 1 and pi^T P = pi^T exactly,
# so the check below is exact in symbolic t1, t2.
t1, t2 = sp.symbols("t1 t2", positive=True)
P_series_1 = sp.eye(n) + Q * t1 + Q**2 * t1**2 / 2
P_series_2 = sp.eye(n) + Q * t2 + Q**2 * t2**2 / 2
M_series = Pim12 * P_series_1.T * Pi * P_series_2 * Pim12
check(
    "M sqrt(pi) = sqrt(pi) exactly",
    zero_matrix(sp.expand(M_series * sqrt_pi - sqrt_pi)),
)
check(
    "sqrt(pi)^T M = sqrt(pi)^T exactly",
    zero_matrix(sp.expand(sqrt_pi.T * M_series - sqrt_pi.T)),
)

# Numerically accurate M for singular values: sigma_0 = 1, sigma_k <= 1,
# and edge reversal (transposition) keeps the singular values.
t1_value, t2_value = sp.Rational(2, 5), sp.Rational(9, 10)
M_numeric = (
    Pim12.evalf(30)
    * expm_numeric(Q.evalf(30) * t1_value).T
    * Pi.evalf(30)
    * expm_numeric(Q.evalf(30) * t2_value)
    * Pim12.evalf(30)
)
import mpmath

mpmath.mp.dps = 30


def singular_values_numeric(matrix):
    mp_matrix = mpmath.matrix([[float(matrix[i, j]) for j in range(n)]
                               for i in range(n)])
    values = mpmath.svd_r(mp_matrix, compute_uv=False)
    return sorted((sp.Float(values[i], 20) for i in range(n)), reverse=True)


sigmas = singular_values_numeric(M_numeric)
check("sigma_0 = 1", abs(sigmas[0] - 1) < sp.Float("1e-12"))
check("all other sigma_k <= 1", all(sv <= 1 for sv in sigmas[1:]))
sigmas_reversed = singular_values_numeric(M_numeric.T)
check(
    "edge reversal keeps the singular values",
    all(abs(a - b) < sp.Float("1e-12") for a, b in zip(sigmas, sigmas_reversed)),
)

# --------------------------------------------------------------------------
# Derivation 3: null variance with cross-mode covariances
# --------------------------------------------------------------------------
print("\n-- 3. Null variance with cross-mode covariances --")

# Right eigenvectors of the two nonstationary modes, exact and complex.
modes = []
for eigenvalue, _multiplicity, vectors in Q.eigenvects():
    if sp.simplify(eigenvalue) != 0:
        modes.append(vectors[0])
check("two nonstationary right eigenvectors", len(modes) == 2)

# Null mean: pi^T r_k = 0 by biorthogonality with the stationary left
# eigenvector, so E[F(k)] = 0 and E[C] = 0 survive non-reversibly.
check(
    "pi^T r_k = 0 for both modes",
    all(is_zero((pi_vec.T * r)[0]) for r in modes),
)

# The Gram matrix in L2(pi) is not the identity: cross-mode covariances.
def gram_entry(k, l):
    return sum(pi_vec[state] * modes[k][state] * modes[l][state]
               for state in range(n))


gram = sp.Matrix(2, 2, lambda k, l: gram_entry(k, l))
check(
    "cross-mode Gram entry G_01 is nonzero",
    not is_zero(gram[0, 1]),
)

# Exact null variance by enumeration over independent stationary sides:
#   C(x, y) = sum_k r_k(x) r_k(y),  weights w_k = 1.
# The conjugate pair keeps C real.
coherence = sp.Matrix(
    n, n,
    lambda x, y: sum(modes[k][x] * modes[k][y] for k in range(2)),
)
check(
    "C(x, y) is real for all state pairs",
    all(is_zero(sp.im(coherence[x, y])) for x in range(n) for y in range(n)),
)

mean_C = sum(
    pi_vec[x] * pi_vec[y] * coherence[x, y] for x in range(n) for y in range(n)
)
check("E[C] = 0 under the saturation null", is_zero(mean_C))

variance_enumerated = sum(
    pi_vec[x] * pi_vec[y] * coherence[x, y] ** 2
    for x in range(n)
    for y in range(n)
)
variance_formula = sum(gram[k, l] ** 2 for k in range(2) for l in range(2))
check(
    "Var(C) = sum_{k,l} G_kl^L G_kl^R exactly",
    is_zero(variance_enumerated - variance_formula),
)

# Restoration of the diagonal rule in the singular basis: for ANY Euclidean
# orthonormal vectors u_k, the L2(pi) Gram of Pi^(-1/2) u_k is the identity.
raw_vectors = [sp.Matrix([1, 2, 0]), sp.Matrix([0, 1, 3]), sp.Matrix([1, 0, 1])]
orthonormal = sp.GramSchmidt(raw_vectors, True)
gram_singular = sp.Matrix(
    n, n,
    lambda k, l: sum(
        pi_vec[state]
        * (orthonormal[k][state] / sp.sqrt(pi_vec[state]))
        * (orthonormal[l][state] / sp.sqrt(pi_vec[state]))
        for state in range(n)
    ),
)
check(
    "L2(pi) Gram of Pi^(-1/2) u_k is the identity",
    zero_matrix(gram_singular - sp.eye(n)),
)

print("\nAll symbolic checks passed.")
