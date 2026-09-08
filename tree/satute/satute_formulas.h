/***************************************************************************
 *   Copyright (C) 2026 by                                                *
 *   IQ-TREE developers                                                    *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 ***************************************************************************/

#ifndef SATUTE_FORMULAS_H
#define SATUTE_FORMULAS_H

#include "satute_types.h"
#include "../phylotree.h"

#include <vector>

namespace satute {
namespace formulas {

/** Selection of substitution-model eigen-modes and the two SatuTe formulas. */
int findZeroEigenIndex(const double *eval, int nstates);
std::vector<int> findNonZeroModes(const double *eval, int nstates, int zero_index);
std::vector<SatuTeFormulaSpec> buildFormulaSpecs(const double *eval, int nstates, int zero_index);
std::vector<double> buildFormulaWeights(
    const double *eval,
    const SatuTeFormulaSpec &formula,
    double effective_branch_length);

/** The two formula labels without model-wide modes, for the non-reversible
 * path where the mode set is per (branch, rate category). */
std::vector<SatuTeFormulaSpec> buildDirectedFormulaSpecs();

/**
 * M(t) = Pi^{1/2} P(t) Pi^{-1/2}, row-major, together with pi. Shared by the
 * directed basis and the directed information fraction.
 */
void buildSymmetrizedOperator(
    ModelSubst *model,
    int nstates,
    double effective_length,
    std::vector<double> &symmetrized,
    std::vector<double> &pi);

/** Build the directed singular basis of M(t) for one effective length. */
void buildDirectedBasis(
    ModelSubst *model,
    int nstates,
    double effective_length,
    SatuTeDirectedBasis &basis);

/**
 * Offsets into the basis of the dominant modes: every sigma within the same
 * relative tolerance of the largest one that the reversible dominant-mode
 * selection uses for near-degenerate eigenvalues.
 */
std::vector<int> dominantModeOffsets(const SatuTeDirectedBasis &basis);

/** Offsets of every retained non-stationary mode. */
std::vector<int> allModeOffsets(const SatuTeDirectedBasis &basis);

/**
 * Formula weights over the selected offsets: all ones for dominant, and
 * sigma_k / sigma_max for the sigma-weighted formula -- the exact analogue of
 * exp((lambda_k - lambda_max) t), to which it reduces reversibly.
 */
std::vector<double> buildSigmaWeights(
    const SatuTeDirectedBasis &basis,
    const std::vector<int> &offsets,
    bool sigma_weighted);

} // namespace formulas
} // namespace satute

#endif
