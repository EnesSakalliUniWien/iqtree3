/***************************************************************************
 *   Copyright (C) 2026 by                                                *
 *   IQ-TREE developers                                                    *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 ***************************************************************************/

#include "satute_formulas.h"
#include "satute_statistics.h"

#include "../phylotree.h"

#include <algorithm>
#include <cmath>
#include <limits>

using namespace std;

namespace satute {
namespace formulas {

int findZeroEigenIndex(const double *eval, int nstates) {
    int zero_index = 0;
    double best_abs = fabs(eval[0]);
    for (int i = 1; i < nstates; i++) {
        double value = fabs(eval[i]);
        if (value < best_abs) {
            best_abs = value;
            zero_index = i;
        }
    }
    return zero_index;
}

vector<int> findNonZeroModes(const double *eval, int nstates, int zero_index) {
    vector<int> modes;
    const double zero_tol = 1e-10;
    for (int i = 0; i < nstates; i++) {
        if (i == zero_index || fabs(eval[i]) <= zero_tol)
            continue;
        modes.push_back(i);
    }
    return modes;
}

static vector<int> findDominantModes(const double *eval, int nstates, int zero_index) {
    vector<int> modes;
    double lambda = -numeric_limits<double>::infinity();
    const double zero_tol = 1e-10;

    for (int i = 0; i < nstates; i++) {
        if (i == zero_index || fabs(eval[i]) <= zero_tol)
            continue;
        if (eval[i] > lambda)
            lambda = eval[i];
    }

    if (!isfinite(lambda))
        return modes;

    double tol = max(1e-8, fabs(lambda) * 1e-6);
    for (int i = 0; i < nstates; i++) {
        if (i == zero_index)
            continue;
        if (fabs(eval[i] - lambda) <= tol)
            modes.push_back(i);
    }
    return modes;
}

vector<SatuTeFormulaSpec> buildFormulaSpecs(const double *eval, int nstates, int zero_index) {
    vector<int> dominant_modes = findDominantModes(eval, nstates, zero_index);
    vector<int> all_modes = findNonZeroModes(eval, nstates, zero_index);
    if (dominant_modes.empty() || all_modes.empty())
        outError("SatuTe cannot identify non-zero substitution-model eigenvalues");

    vector<SatuTeFormulaSpec> formulas;
    SatuTeFormulaSpec dominant;
    dominant.name = "dominant";
    dominant.modes = dominant_modes;
    dominant.eigenvalue_weighted = false;
    formulas.push_back(dominant);

    SatuTeFormulaSpec weighted;
    weighted.name = "eigenvalue_weighted";
    weighted.modes = all_modes;
    weighted.eigenvalue_weighted = true;
    formulas.push_back(weighted);
    return formulas;
}

vector<double> buildFormulaWeights(
    const double *eval,
    const SatuTeFormulaSpec &formula,
    double effective_branch_length)
{
    vector<double> weights(formula.modes.size(), 1.0);
    if (!formula.eigenvalue_weighted)
        return weights;

    double dominant = -numeric_limits<double>::infinity();
    for (size_t i = 0; i < formula.modes.size(); i++)
        dominant = max(dominant, eval[formula.modes[i]]);

    for (size_t i = 0; i < formula.modes.size(); i++)
        weights[i] = exp((eval[formula.modes[i]] - dominant) * effective_branch_length);
    return weights;
}

vector<SatuTeFormulaSpec> buildDirectedFormulaSpecs() {
    vector<SatuTeFormulaSpec> formulas;
    SatuTeFormulaSpec dominant;
    dominant.name = "dominant";
    dominant.eigenvalue_weighted = false;
    formulas.push_back(dominant);

    SatuTeFormulaSpec weighted;
    weighted.name = "eigenvalue_weighted";
    weighted.eigenvalue_weighted = true;
    formulas.push_back(weighted);
    return formulas;
}

/**
 * M(t) = Pi^{1/2} P(t) Pi^{-1/2}, row-major, together with pi. The transition
 * matrix comes from computeTransMatrix(), which dispatches to the
 * non-reversible matrix exponential -- the same operator the likelihood
 * kernel uses, so the statistic tests exactly the fitted transition law.
 */
void buildSymmetrizedOperator(
    ModelSubst *model,
    int nstates,
    double effective_length,
    vector<double> &symmetrized,
    vector<double> &pi)
{
    pi.assign(nstates, 0.0);
    model->getStateFrequency(&pi[0]);
    for (int state = 0; state < nstates; state++) {
        if (!(pi[state] > 0.0))
            outError("SatuTe requires strictly positive stationary frequencies for non-reversible models");
    }

    vector<double> transition((size_t)nstates * nstates, 0.0);
    model->computeTransMatrix(effective_length, &transition[0]);

    symmetrized.assign((size_t)nstates * nstates, 0.0);
    for (int x = 0; x < nstates; x++) {
        for (int y = 0; y < nstates; y++) {
            symmetrized[(size_t)x * nstates + y] =
                sqrt(pi[x] / pi[y]) * transition[(size_t)x * nstates + y];
        }
    }
}

void buildDirectedBasis(
    ModelSubst *model,
    int nstates,
    double effective_length,
    SatuTeDirectedBasis &basis)
{
    basis.nstates = nstates;
    basis.sigma.clear();
    basis.left_projection.clear();
    basis.right_projection.clear();

    vector<double> symmetrized;
    vector<double> pi;
    buildSymmetrizedOperator(model, nstates, effective_length, symmetrized, pi);

    // Singular value decomposition through the symmetric eigenproblem of
    // M^T M = V diag(sigma^2) V^T, then u_k = M v_k / sigma_k. The matrices
    // are tiny (4x4 for DNA), so the shared Jacobi solver is adequate.
    vector<double> mtm((size_t)nstates * nstates, 0.0);
    for (int i = 0; i < nstates; i++) {
        for (int j = 0; j < nstates; j++) {
            double value = 0.0;
            for (int x = 0; x < nstates; x++) {
                value += symmetrized[(size_t)x * nstates + i] *
                    symmetrized[(size_t)x * nstates + j];
            }
            mtm[(size_t)i * nstates + j] = value;
        }
    }
    vector<double> right_vectors;   // columns of V
    jacobiEigenSymmetric(mtm, nstates, right_vectors);

    struct DirectedMode {
        double sigma;
        vector<double> u;
        vector<double> v;
    };
    vector<DirectedMode> modes;
    int stationary_index = -1;
    double stationary_overlap = -1.0;

    vector<double> sqrt_pi(nstates, 0.0);
    for (int state = 0; state < nstates; state++)
        sqrt_pi[state] = sqrt(pi[state]);

    for (int k = 0; k < nstates; k++) {
        DirectedMode mode;
        mode.sigma = sqrt(max(0.0, mtm[(size_t)k * nstates + k]));
        mode.v.assign(nstates, 0.0);
        for (int state = 0; state < nstates; state++)
            mode.v[state] = right_vectors[(size_t)state * nstates + k];

        mode.u.assign(nstates, 0.0);
        if (mode.sigma > 1e-12) {
            for (int x = 0; x < nstates; x++) {
                double value = 0.0;
                for (int y = 0; y < nstates; y++)
                    value += symmetrized[(size_t)x * nstates + y] * mode.v[y];
                mode.u[x] = value / mode.sigma;
            }
        }

        double overlap = 0.0;
        for (int state = 0; state < nstates; state++)
            overlap += mode.v[state] * sqrt_pi[state];
        overlap = fabs(overlap);

        modes.push_back(mode);
        if (overlap > stationary_overlap) {
            stationary_overlap = overlap;
            stationary_index = (int)modes.size() - 1;
        }
    }

    // Drop the stationary triplet (sigma = 1, u = v = sqrt(pi)); keep the
    // rest sorted by descending canonical correlation, and drop numerically
    // vanished modes whose left vectors are undefined -- their weight and
    // contribution are zero at working precision anyway.
    vector<int> order;
    for (int k = 0; k < (int)modes.size(); k++) {
        if (k == stationary_index || modes[k].sigma <= 1e-12)
            continue;
        order.push_back(k);
    }
    stable_sort(order.begin(), order.end(), [&modes](int a, int b) {
        return modes[a].sigma > modes[b].sigma;
    });

    for (size_t rank = 0; rank < order.size(); rank++) {
        const DirectedMode &mode = modes[order[rank]];
        basis.sigma.push_back(mode.sigma);
        for (int state = 0; state < nstates; state++)
            basis.left_projection.push_back(mode.u[state] / sqrt_pi[state]);
        for (int state = 0; state < nstates; state++)
            basis.right_projection.push_back(mode.v[state] / sqrt_pi[state]);
    }
}

vector<int> dominantModeOffsets(const SatuTeDirectedBasis &basis) {
    vector<int> offsets;
    if (basis.sigma.empty())
        return offsets;
    double largest = basis.sigma[0];
    double tolerance = max(1e-8, largest * 1e-6);
    for (size_t k = 0; k < basis.sigma.size(); k++) {
        if (fabs(basis.sigma[k] - largest) <= tolerance)
            offsets.push_back((int)k);
    }
    return offsets;
}

vector<int> allModeOffsets(const SatuTeDirectedBasis &basis) {
    vector<int> offsets;
    for (size_t k = 0; k < basis.sigma.size(); k++)
        offsets.push_back((int)k);
    return offsets;
}

vector<double> buildSigmaWeights(
    const SatuTeDirectedBasis &basis,
    const vector<int> &offsets,
    bool sigma_weighted)
{
    vector<double> weights(offsets.size(), 1.0);
    if (!sigma_weighted || offsets.empty())
        return weights;

    double largest = 0.0;
    for (size_t i = 0; i < offsets.size(); i++)
        largest = max(largest, basis.sigma[offsets[i]]);
    if (!(largest > 0.0))
        return weights;
    for (size_t i = 0; i < offsets.size(); i++)
        weights[i] = basis.sigma[offsets[i]] / largest;
    return weights;
}

} // namespace formulas
} // namespace satute
