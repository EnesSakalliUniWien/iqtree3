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

} // namespace formulas
} // namespace satute

#endif
