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

#include <memory>
#include <vector>

namespace satute {

enum class StatisticKind {
    Category,
    Mixture
};

struct StatisticRequest {
    Node *left_node;
    Node *left_dad;
    Node *right_node;
    Node *right_dad;
    Alignment *aln;
    ModelSubst *model;
    int nstates;
    const double *eigenvalues;
    const SatuTeFormulaSpec *formula;
    const std::vector<double> *mode_weights;
    const IntVector *pattern_categories;
    int required_pattern_category;
    double rate_multiplier;
    const std::vector<int> *modes;
    const std::vector<SatuTeRateCategory> *rate_categories;
    double branch_length;
    double log_weight_shift;
    const SatuTeBranchResult *base;
};

class StatisticStrategy {
public:
    virtual ~StatisticStrategy() {}
    virtual SatuTeBranchResult compute(const StatisticRequest &request) const = 0;
};

class StatisticStrategyFactory {
public:
    static std::unique_ptr<StatisticStrategy> create(StatisticKind kind);
};

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

SatuTeBranchResult computeMixtureStatistic(
    Node *left_node,
    Node *left_dad,
    Node *right_node,
    Node *right_dad,
    Alignment *aln,
    ModelSubst *model,
    int nstates,
    const double *eval,
    const std::vector<int> &modes,
    const std::vector<SatuTeRateCategory> &categories,
    double branch_length,
    double log_weight_shift,
    const SatuTeBranchResult &base);

SatuTeBranchResult poolCategoryResults(
    const std::vector<SatuTeBranchResult> &category_results,
    const SatuTeBranchResult &base);

} // namespace satute

#endif
