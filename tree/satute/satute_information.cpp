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
#include "satute_formulas.h"

#include <cmath>
#include <functional>
#include <limits>

using namespace std;

namespace satute {
namespace information {

/** Rate-proportion-weighted pooling shared by both information fractions. */
static double poolInformationOverCategories(
    const vector<SatuTeRateCategory> &categories,
    double branch_length,
    const std::function<double(double)> &category_information)
{
    double weighted_information = 0.0;
    double total_proportion = 0.0;
    for (size_t i = 0; i < categories.size(); i++) {
        const SatuTeRateCategory &category = categories[i];
        if (!(category.proportion > 0.0) || !isfinite(category.proportion))
            continue;
        double information = category_information(branch_length * category.rate);
        if (!isfinite(information))
            return numeric_limits<double>::quiet_NaN();
        weighted_information += category.proportion * information;
        total_proportion += category.proportion;
    }
    if (!(total_proportion > 0.0))
        return numeric_limits<double>::quiet_NaN();
    return weighted_information / total_proportion;
}

/** Clamp tiny numerical excursions outside [0, 1]. */
static double clampInformation(double information) {
    if (information < 0.0 && information > -1e-12)
        information = 0.0;
    if (information > 1.0 && information < 1.0 + 1e-12)
        information = 1.0;
    return information;
}

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
    return clampInformation(information);
}

double computePooledInformationFraction(
    const double *eval,
    const vector<int> &modes,
    const vector<SatuTeRateCategory> &categories,
    double branch_length)
{
    return poolInformationOverCategories(
        categories,
        branch_length,
        [eval, &modes](double effective_length) {
            return computeInformationFraction(eval, modes, effective_length);
        });
}

double computeDirectedInformationFraction(
    ModelSubst *model,
    int nstates,
    double effective_length)
{
    if (nstates < 2 || effective_length < 0.0 || !isfinite(effective_length))
        return numeric_limits<double>::quiet_NaN();

    vector<double> symmetrized;
    vector<double> pi;
    formulas::buildSymmetrizedOperator(
        model, nstates, effective_length, symmetrized, pi);

    double energy = 0.0;
    for (int x = 0; x < nstates; x++) {
        for (int y = 0; y < nstates; y++) {
            double stationary_component = sqrt(pi[x] * pi[y]);
            double difference =
                symmetrized[(size_t)x * nstates + y] - stationary_component;
            energy += difference * difference;
        }
    }
    return clampInformation(energy / (nstates - 1));
}

double computePooledDirectedInformationFraction(
    ModelSubst *model,
    int nstates,
    const vector<SatuTeRateCategory> &categories,
    double branch_length)
{
    return poolInformationOverCategories(
        categories,
        branch_length,
        [model, nstates](double effective_length) {
            return computeDirectedInformationFraction(model, nstates, effective_length);
        });
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
