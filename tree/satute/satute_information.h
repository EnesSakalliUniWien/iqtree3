/***************************************************************************
 *   Copyright (C) 2026 by                                                *
 *   IQ-TREE developers                                                    *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 ***************************************************************************/

#ifndef SATUTE_INFORMATION_H
#define SATUTE_INFORMATION_H

#include "satute_types.h"

#include <vector>

class ModelSubst;

namespace satute {
namespace information {

/** Retained-information fraction per branch and its saturation-index scaling. */
double computeInformationFraction(
    const double *eval,
    const std::vector<int> &modes,
    double effective_length);
double computePooledInformationFraction(
    const double *eval,
    const std::vector<int> &modes,
    const std::vector<SatuTeRateCategory> &categories,
    double branch_length);

/**
 * Non-reversible information fraction: the normalized squared Frobenius
 * energy ||M(t) - Pi0||_F^2 / K with K = nstates - 1, where
 * M(t) = Pi^{1/2} P(t) Pi^{-1/2} and Pi0 = sqrt(pi) sqrt(pi)^T. The
 * eigenvalue-sum formula is invalid without an orthonormal eigenbasis, but
 * this operator definition reduces to it reversibly and remains monotone
 * non-reversibly, because only the symmetric part of the symmetrized
 * generator dissipates; the antisymmetric part rotates without dissipating.
 */
double computeDirectedInformationFraction(
    ModelSubst *model,
    int nstates,
    double effective_length);
double computePooledDirectedInformationFraction(
    ModelSubst *model,
    int nstates,
    const std::vector<SatuTeRateCategory> &categories,
    double branch_length);

void setSaturationScale(SatuTeBranchResult &result, double information_fraction);

} // namespace information
} // namespace satute

#endif
