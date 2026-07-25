/***************************************************************************
 *   Copyright (C) 2026 by                                                *
 *   IQ-TREE developers                                                    *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 ***************************************************************************/

#include "satute_results.h"

#include <algorithm>
#include <limits>

using namespace std;

namespace satute {
namespace results {

SatuTeBranchResult makeBaseResult(
    int branch_id,
    int left_taxa,
    int right_taxa,
    double branch_length,
    double effective_length,
    double alpha,
    const string &formula,
    const string &rate_category,
    const string &rate_multiplier,
    const string &modes,
    const string &eigenvalues,
    const string &weights,
    const string &label,
    const string &split)
{
    SatuTeBranchResult result;
    result.id = branch_id;
    result.left_taxa = left_taxa;
    result.right_taxa = right_taxa;
    result.valid_sites = 0;
    result.skipped_sites = 0;
    result.length = branch_length;
    result.effective_length = effective_length;
    result.information_fraction = numeric_limits<double>::quiet_NaN();
    result.saturation_index = numeric_limits<double>::quiet_NaN();
    result.coherence = numeric_limits<double>::quiet_NaN();
    result.variance = numeric_limits<double>::quiet_NaN();
    result.se = numeric_limits<double>::quiet_NaN();
    result.z_score = numeric_limits<double>::quiet_NaN();
    result.p_value = numeric_limits<double>::quiet_NaN();
    result.fdr_by = numeric_limits<double>::quiet_NaN();
    result.alpha = alpha;
    result.alpha_adjusted = alpha / max(1.0, (double)left_taxa * (double)right_taxa);
    result.formula = formula;
    result.rate_category = rate_category;
    result.rate_multiplier = rate_multiplier;
    result.rate_category_sites = 0;
    result.modes = modes;
    result.eigenvalues = eigenvalues;
    result.weights = weights;
    result.decision = "undefined";
    result.decision_bonferroni = "undefined";
    result.decision_fdr = "not_tested";
    result.label = label;
    result.split = split;
    return result;
}

} // namespace results
} // namespace satute
