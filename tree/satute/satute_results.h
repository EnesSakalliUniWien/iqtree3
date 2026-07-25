/***************************************************************************
 *   Copyright (C) 2026 by                                                *
 *   IQ-TREE developers                                                    *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 ***************************************************************************/

#ifndef SATUTE_RESULTS_H
#define SATUTE_RESULTS_H

#include "satute_types.h"

#include <string>

namespace satute {
namespace results {

/** Construct a fully-initialised result row before any statistic is computed. */
SatuTeBranchResult makeBaseResult(
    int branch_id,
    int left_taxa,
    int right_taxa,
    double branch_length,
    double effective_length,
    double alpha,
    const std::string &formula,
    const std::string &rate_category,
    const std::string &rate_multiplier,
    const std::string &modes,
    const std::string &eigenvalues,
    const std::string &weights,
    const std::string &label,
    const std::string &split);

} // namespace results
} // namespace satute

#endif
