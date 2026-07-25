/***************************************************************************
 *   Copyright (C) 2026 by                                                *
 *   IQ-TREE developers                                                    *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 ***************************************************************************/

#ifndef SATUTE_RATES_H
#define SATUTE_RATES_H

#include "satute_types.h"
#include "../phylotree.h"

#include <vector>

namespace satute {
namespace rates {

/** Map IQ-TREE's fitted rate model onto the per-category rows SatuTe reports. */
std::vector<SatuTeRateCategory> buildRateCategories(
    RateHeterogeneity *site_rate,
    DoubleVector &pattern_rates,
    IntVector &pattern_cat);

} // namespace rates
} // namespace satute

#endif
