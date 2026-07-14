/***************************************************************************
 *   Copyright (C) 2026 by                                                *
 *   IQ-TREE developers                                                    *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 ***************************************************************************/

#ifndef SATUTE_STATISTICS_H
#define SATUTE_STATISTICS_H

#include "satute_types.h"
#include "../phylotree.h"

#include <vector>

namespace satute {

SatuTeBranchResult computeCategoryStatistic(
    Node *left_node,
    Node *left_dad,
    Node *right_node,
    Node *right_dad,
    Alignment *aln,
    ModelSubst *model,
    int nstates,
    const double *eval,
    const SatuTeFormulaSpec &formula,
    const std::vector<double> &mode_weights,
    const IntVector *pattern_cat,
    int required_pattern_cat,
    double rate_multiplier,
    const SatuTeBranchResult &base);

SatuTeBranchResult poolCategoryResults(
    const std::vector<SatuTeBranchResult> &category_results,
    const SatuTeBranchResult &base);

/**
 * Apply Benjamini-Yekutieli FDR adjustment independently to the pooled
 * dominant and eigenvalue_weighted branch p-values. Rate-category rows are
 * diagnostic and are not members of either FDR family.
 */
void applySeparateFormulaFdr(
    std::vector<SatuTeBranchResult> &results,
    double fdr_level);

} // namespace satute

#endif
