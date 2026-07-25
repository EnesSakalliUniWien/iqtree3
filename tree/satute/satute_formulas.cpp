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

} // namespace formulas
} // namespace satute
