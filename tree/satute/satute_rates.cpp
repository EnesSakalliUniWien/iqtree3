/***************************************************************************
 *   Copyright (C) 2026 by                                                *
 *   IQ-TREE developers                                                    *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 ***************************************************************************/

#include "satute_rates.h"

#include <cmath>

using namespace std;

namespace satute {
namespace rates {

vector<SatuTeRateCategory> buildRateCategories(
    RateHeterogeneity *site_rate,
    DoubleVector &pattern_rates,
    IntVector &pattern_cat)
{
    vector<SatuTeRateCategory> categories;
    bool single_rate =
        site_rate->getNRate() == 1 &&
        fabs(site_rate->getRate(0) - 1.0) <= 1e-8 &&
        fabs(site_rate->getProp(0) - 1.0) <= 1e-8 &&
        fabs(site_rate->getPInvar()) <= 1e-12;

    if (single_rate) {
        SatuTeRateCategory category;
        category.label = "pooled";
        category.pattern_cat = -1;
        category.rate = 1.0;
        category.proportion = 1.0;
        categories.push_back(category);
        return categories;
    }

    // Keep rate-category assignment identical to IQ-TREE's empirical-Bayes
    // site-rate path: best posterior category per pattern, with IQ-TREE's
    // fitted category multipliers from getRate().
    int observed_categories = site_rate->computePatternRates(pattern_rates, pattern_cat);
    if (pattern_cat.empty()) {
        if (site_rate->getPInvar() > 0.0) {
            outError("SatuTe rate-category analysis requires IQ-TREE pattern-rate assignments. "
                     "IQ-TREE does not emit per-site category rows for a pure +I model; use +I with +G/+R or a discrete rate model.");
        }
        outError("SatuTe could not obtain IQ-TREE pattern-rate assignments for this rate model");
    }

    if (site_rate->getPInvar() > 0.0) {
        SatuTeRateCategory invariant;
        invariant.label = "0";
        invariant.pattern_cat = 0;
        invariant.rate = 0.0;
        invariant.proportion = site_rate->getPInvar();
        categories.push_back(invariant);
        for (int c = 0; c < site_rate->getNRate(); c++) {
            SatuTeRateCategory category;
            category.label = convertIntToString(c + 1);
            category.pattern_cat = c + 1;
            category.rate = site_rate->getRate(c);
            category.proportion = site_rate->getProp(c);
            categories.push_back(category);
        }
    } else {
        for (int c = 0; c < site_rate->getNRate(); c++) {
            SatuTeRateCategory category;
            category.label = convertIntToString(c + 1);
            category.pattern_cat = c;
            category.rate = site_rate->getRate(c);
            category.proportion = site_rate->getProp(c);
            categories.push_back(category);
        }
    }

    if (observed_categories != (int)categories.size()) {
        outWarning("SatuTe rate-category count differs from IQ-TREE pattern-rate assignment count");
    }
    return categories;
}

} // namespace rates
} // namespace satute
