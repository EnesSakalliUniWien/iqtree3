/***************************************************************************
 *   Copyright (C) 2026 by                                                *
 *   IQ-TREE developers                                                    *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 ***************************************************************************/

#include "satute_information.h"

#include <cmath>
#include <limits>

using namespace std;

namespace satute {
namespace information {

double computeInformationFraction(
    const double *eval,
    const vector<int> &modes,
    double effective_length)
{
    if (modes.empty() || effective_length < 0.0 || !isfinite(effective_length))
        return numeric_limits<double>::quiet_NaN();

    double information = 0.0;
    for (size_t i = 0; i < modes.size(); i++)
        information += exp(2.0 * eval[modes[i]] * effective_length);
    information /= modes.size();
    if (information < 0.0 && information > -1e-12)
        information = 0.0;
    if (information > 1.0 && information < 1.0 + 1e-12)
        information = 1.0;
    return information;
}

double computePooledInformationFraction(
    const double *eval,
    const vector<int> &modes,
    const vector<SatuTeRateCategory> &categories,
    double branch_length)
{
    double weighted_information = 0.0;
    double total_proportion = 0.0;
    for (size_t i = 0; i < categories.size(); i++) {
        const SatuTeRateCategory &category = categories[i];
        if (!(category.proportion > 0.0) || !isfinite(category.proportion))
            continue;
        double category_information = computeInformationFraction(
            eval,
            modes,
            branch_length * category.rate);
        if (!isfinite(category_information))
            return numeric_limits<double>::quiet_NaN();
        weighted_information += category.proportion * category_information;
        total_proportion += category.proportion;
    }
    if (!(total_proportion > 0.0))
        return numeric_limits<double>::quiet_NaN();
    return weighted_information / total_proportion;
}

void setSaturationScale(SatuTeBranchResult &result, double information_fraction) {
    result.information_fraction = information_fraction;
    result.saturation_index = isfinite(information_fraction)
        ? 1.0 - information_fraction
        : numeric_limits<double>::quiet_NaN();
    if (result.saturation_index < 0.0 && result.saturation_index > -1e-12)
        result.saturation_index = 0.0;
    if (result.saturation_index > 1.0 && result.saturation_index < 1.0 + 1e-12)
        result.saturation_index = 1.0;
}

} // namespace information
} // namespace satute
