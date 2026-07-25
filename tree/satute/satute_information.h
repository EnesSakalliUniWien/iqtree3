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
void setSaturationScale(SatuTeBranchResult &result, double information_fraction);

} // namespace information
} // namespace satute

#endif
