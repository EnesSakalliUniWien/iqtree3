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

/**
 * Directional partial likelihoods: for every (node, neighbour) pair, the
 * message L(node \ neighbour) at every alignment pattern, for ONE rate
 * category.
 *
 * Why this exists: the per-branch recursion recomputes an entire subtree for
 * each of the 2n-3 branches, so the statistic was O(n^2) in taxa. Each branch
 * needs only the two messages that meet on it, and the whole set of 2(2n-3)
 * messages is produced by two traversals per pattern. Filling this once per
 * rate category makes the statistic O(n) in taxa.
 *
 * Keyed by node id and the neighbour's index in that node's neighbour list --
 * dense, and independent of any rooting.
 */
struct SatuTeDirectionalCache {
    std::vector<double> data;
    int max_degree = 0;
    int nstates = 0;
    size_t nptn = 0;
    size_t dir_stride = 0;   // nptn * nstates
    bool ready = false;

    inline const double *message(int node_id, int neighbor_index, size_t ptn) const {
        return &data[((size_t)node_id * (size_t)max_degree + (size_t)neighbor_index)
                     * dir_stride + ptn * (size_t)nstates];
    }
    inline double *message(int node_id, int neighbor_index, size_t ptn) {
        return &data[((size_t)node_id * (size_t)max_degree + (size_t)neighbor_index)
                     * dir_stride + ptn * (size_t)nstates];
    }
};

/** Bytes a cache of this shape needs, for budgeting before allocating. */
size_t directionalCacheBytes(
    size_t node_bound, int max_degree, size_t nptn, int nstates);

/**
 * Fill every directional message for one rate category. Returns false when the
 * tree has no interior node to traverse from, in which case the caller keeps
 * using the per-branch recursion.
 */
bool buildDirectionalCache(
    Node *root,
    Alignment *aln,
    ModelSubst *model,
    int nstates,
    double rate_multiplier,
    int num_threads,
    size_t node_bound,
    int max_degree,
    SatuTeDirectionalCache &cache);

/**
 * @param out_pattern_coherence  Optional. When non-null, resized to the
 * pattern count and filled with each pattern's centred, weighted contribution
 * to satC:  sum_i w_i (L_i - mean_L_i)(R_i - mean_R_i).  These sum, frequency
 * weighted and divided by valid_sites, to exactly the returned result's
 * `coherence` -- the invariant the .sat.site writer is checked against.
 * Patterns outside this call's rate category (or skipped) are left at 0, so
 * per-category arrays tile the alignment and can simply be summed for pooled.
 */
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
    int num_threads,
    bool enhanced_weights,
    const SatuTeBranchResult &base,
    const SatuTeDirectionalCache *directional_cache = nullptr,
    std::vector<double> *out_pattern_coherence = nullptr);

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
