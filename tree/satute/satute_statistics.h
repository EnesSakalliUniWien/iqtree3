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
    std::vector<double> log_scales;
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
    inline double &logScale(int node_id, int neighbor_index, size_t ptn) {
        return log_scales[((size_t)node_id * max_degree + neighbor_index) * nptn + ptn];
    }
    inline double logScale(int node_id, int neighbor_index, size_t ptn) const {
        return log_scales[((size_t)node_id * max_degree + neighbor_index) * nptn + ptn];
    }
};

/** Streaming centred mixture score, with a common persistence scale. */
struct SatuTeMixtureAccumulator {
    double proportion = 0.0;
    double weight_scale = 0.0;
    std::vector<double> log_null;
    std::vector<double> score;
    std::vector<char> invalid;
};

/** Bytes a cache of this shape needs, for budgeting before allocating. */
size_t directionalCacheBytes(
    size_t node_bound, int max_degree, size_t nptn, int nstates,
    bool track_scales = false);

/**
 * Symmetric eigenvalue decomposition by cyclic Jacobi rotations: a, IN/OUT,
 * becomes diagonal; eigenvectors, OUT, has the eigenvectors as columns.
 * Adequate for the small matrices here (n <= ~19 for protein, typically 1-3
 * for DNA): a handful of sweeps converges to machine precision at this size,
 * and no external linear-algebra dependency is needed for it. Shared with the
 * non-reversible directed-basis builder, which diagonalizes M^T M for its
 * singular value decomposition.
 */
void jacobiEigenSymmetric(
    std::vector<double> &a, int n, std::vector<double> &eigenvectors);

/**
 * Directed (rooted, non-reversible) evaluation context for
 * computeCategoryStatistic. When supplied:
 * - side partial likelihoods are recomputed direction-aware on the rooted
 *   tree: edges toward the root transport with the transposed transition
 *   matrix and the root leaf contributes the stationary distribution, so the
 *   rootward side factor is the joint P(rootward data, X_u = x) and is NOT
 *   multiplied by the stationary frequencies again;
 * - side factors project onto the left/right singular bases of the
 *   symmetrized joint operator M = Pi^{1/2} P(t) Pi^{-1/2} instead of the
 *   shared reversible eigenvector basis.
 */
struct SatuTeDirectedContext {
    const double *left_projection;   // nmodes x nstates: u_k(x)/sqrt(pi_x)
    const double *right_projection;  // nmodes x nstates: v_k(y)/sqrt(pi_y)
    const std::vector<int> *parent_of;  // node id -> parent node id; -1 at root
    int root_node_id;                   // id of the root leaf carrying pi
};

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
    SatuTeDirectionalCache &cache,
    bool track_scales = false);

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
    std::vector<double> *out_pattern_coherence = nullptr,
    const SatuTeDirectedContext *directed = nullptr,
    SatuTeMixtureAccumulator *mixture = nullptr);

SatuTeBranchResult computeMixtureStatistic(
    const SatuTeMixtureAccumulator &mixture,
    Alignment *aln,
    const SatuTeBranchResult &base,
    std::vector<double> *out_pattern_coherence = nullptr);

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
