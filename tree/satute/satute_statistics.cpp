/***************************************************************************
 *   Copyright (C) 2026 by                                                *
 *   IQ-TREE developers                                                    *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 ***************************************************************************/

#include "satute_statistics.h"
#include "satute_format.h"
#include "../../gsl/mygsl.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <limits>
#include <map>
#include <utility>
#include <vector>

using namespace std;

namespace satute {

/**
 * Transition matrices keyed by the id of the node the branch descends into.
 * A directional traversal reaches each edge once, and the matrix depends only
 * on branch length, so that id is a stable dense key. A map<Neighbor*,...> here
 * instead cost a red-black-tree walk per node visit per pattern, which was
 * latency bound and throttled the parallel loop.
 */
static void ensureTransitionMatrix(
    Neighbor *branch,
    ModelSubst *model,
    int nstates,
    double rate_multiplier,
    vector<double> &transition_cache)
{
    double *eigenvectors = model->getEigenvectors();
    double *inv_eigenvectors = model->getInverseEigenvectors();
    double *eval = model->getEigenvalues();
    const int matrix_size = nstates * nstates;
    double *transition = &transition_cache[(size_t)branch->node->id * matrix_size];
    double length = branch->length * rate_multiplier;

    if (length == 0.0) {
        for (int i = 0; i < nstates; i++)
            for (int j = 0; j < nstates; j++)
                transition[i * nstates + j] = (i == j) ? 1.0 : 0.0;
        return;
    }

    for (int parent_state = 0; parent_state < nstates; parent_state++) {
        for (int child_state = 0; child_state < nstates; child_state++) {
            double value = 0.0;
            for (int eigen = 0; eigen < nstates; eigen++) {
                value += eigenvectors[parent_state * nstates + eigen] *
                    inv_eigenvectors[eigen * nstates + child_state] *
                    exp(eval[eigen] * length);
            }
            transition[parent_state * nstates + child_state] = value;
        }
    }
}

/**
 * Populate the transition-matrix cache for every branch the partial-likelihood
 * recursion below will visit. The traversal mirrors computeOrdinaryPartial()
 * exactly. Matrices depend only on the branch and the rate multiplier, never on
 * the alignment pattern, so filling them up front turns the cache into
 * read-only shared state that the parallel pattern loop can share safely.
 */
static void prepareTransitionMatrices(
    Node *node,
    Node *dad,
    ModelSubst *model,
    int nstates,
    double rate_multiplier,
    vector<double> &transition_cache)
{
    if (node->isLeaf())
        return;

    FOR_NEIGHBOR_IT(node, dad, it) {
        Neighbor *child_branch = *it;
        ensureTransitionMatrix(child_branch, model, nstates, rate_multiplier, transition_cache);
        prepareTransitionMatrices(
            child_branch->node,
            node,
            model,
            nstates,
            rate_multiplier,
            transition_cache);
    }
}

static void computeOrdinaryPartial(
    Node *node,
    Node *dad,
    size_t ptn,
    Alignment *aln,
    ModelSubst *model,
    int nstates,
    const vector<double> &transition_cache,
    double *scratch,
    double *partial,
    double *log_scale = nullptr)
{
    if (node->isLeaf()) {
        // IUPAC ambiguity codes occupy states in [nstates, STATE_UNKNOWN) and
        // must be expanded to the subset of states they stand for. Let the model
        // do it, exactly as PhyloTree::computeTipPartialLikelihood() does, so
        // that model-specific handling (ModelDNA, ModelProtein, ModelDNAError)
        // stays in effect. Note the out-of-range fallback must be STATE_UNKNOWN
        // rather than nstates, which ModelDNA would decode as a real ambiguity
        // code.
        StateType state = (node->id >= 0 && node->id < (int)aln->getNSeq())
            ? aln->at(ptn)[node->id]
            : aln->STATE_UNKNOWN;
        model->computeTipLikelihood(state, partial);
        return;
    }

    for (int state = 0; state < nstates; state++)
        partial[state] = 1.0;

    // This level's child buffer is carved from the caller's scratch block; the
    // recursion hands the remainder down. Allocating a vector per node visit
    // here instead put a malloc/free pair in the innermost loop, which under
    // OpenMP became allocator contention and capped the speedup.
    double *child_partial = scratch;
    FOR_NEIGHBOR_IT(node, dad, it) {
        Neighbor *child_branch = *it;
        computeOrdinaryPartial(
            child_branch->node,
            node,
            ptn,
            aln,
            model,
            nstates,
            transition_cache,
            scratch + nstates,
            child_partial,
            log_scale);

        // Filled by prepareTransitionMatrices(), which walks the same edges.
        const double *transition =
            &transition_cache[(size_t)child_branch->node->id * nstates * nstates];
        for (int parent_state = 0; parent_state < nstates; parent_state++) {
            double contribution = 0.0;
            for (int child_state = 0; child_state < nstates; child_state++) {
                contribution += transition[parent_state * nstates + child_state] *
                    child_partial[child_state];
            }
            partial[parent_state] *= contribution;
        }
    }

    // Mixture responsibilities need the likelihood scale that posteriors cancel.
    double partial_sum = 0.0;
    for (int state = 0; state < nstates; state++)
        partial_sum += partial[state];
    if (partial_sum > 0.0 && isfinite(partial_sum)) {
        if (log_scale != nullptr)
            *log_scale += log(partial_sum);
        for (int state = 0; state < nstates; state++)
            partial[state] /= partial_sum;
    }
}

/**
 * Direction-aware transition matrices for the rooted non-reversible path.
 * Each edge keeps ONE matrix, computed for the edge's TRUE direction (parent
 * to child under the tree rooting), keyed by the true child's node id. The
 * traversal mirrors computeDirectedPartial() below. computeTransMatrix()
 * dispatches to the non-reversible matrix exponential for these models.
 */
static void prepareDirectedTransitionMatrices(
    Node *node,
    Node *dad,
    ModelSubst *model,
    int nstates,
    double rate_multiplier,
    const vector<int> &parent_of,
    vector<double> &transition_cache,
    vector<char> &transition_filled)
{
    if (node->isLeaf())
        return;

    FOR_NEIGHBOR_IT(node, dad, it) {
        Neighbor *branch = *it;
        Node *neighbor = branch->node;
        Node *true_child = (parent_of[neighbor->id] == node->id) ? neighbor : node;
        if (!transition_filled[true_child->id]) {
            model->computeTransMatrix(
                branch->length * rate_multiplier,
                &transition_cache[(size_t)true_child->id * nstates * nstates]);
            transition_filled[true_child->id] = 1;
        }
        prepareDirectedTransitionMatrices(
            neighbor, node, model, nstates, rate_multiplier, parent_of,
            transition_cache, transition_filled);
    }
}

/**
 * Direction-aware side vector for the rooted non-reversible path: the product
 * of the messages into `node` from every neighbour except `dad`.
 *
 * Downward edges (neighbour below node) transport a below-likelihood with the
 * parent-to-child matrix contracted over child states, exactly as the
 * reversible recursion does. Upward edges (neighbour above node) transport
 * with the SAME matrix contracted over the PARENT states -- the transpose
 * direction -- and the root leaf contributes the stationary distribution
 * instead of an all-ones tip vector. The result on the rootward side of a
 * focal branch is therefore the joint P(rootward data, X = x), with the
 * stationary prior already carried down from the root.
 */
static void computeDirectedPartial(
    Node *node,
    Node *dad,
    size_t ptn,
    Alignment *aln,
    ModelSubst *model,
    int nstates,
    const vector<int> &parent_of,
    int root_node_id,
    const double *state_freq,
    const vector<double> &transition_cache,
    double *scratch,
    double *partial,
    double *log_scale = nullptr)
{
    if (node->isLeaf()) {
        if (node->id == root_node_id) {
            for (int state = 0; state < nstates; state++)
                partial[state] = state_freq[state];
            return;
        }
        StateType state = (node->id >= 0 && node->id < (int)aln->getNSeq())
            ? aln->at(ptn)[node->id]
            : aln->STATE_UNKNOWN;
        model->computeTipLikelihood(state, partial);
        return;
    }

    for (int state = 0; state < nstates; state++)
        partial[state] = 1.0;

    double *neighbor_partial = scratch;
    FOR_NEIGHBOR_IT(node, dad, it) {
        Neighbor *branch = *it;
        Node *neighbor = branch->node;
        computeDirectedPartial(
            neighbor,
            node,
            ptn,
            aln,
            model,
            nstates,
            parent_of,
            root_node_id,
            state_freq,
            transition_cache,
            scratch + nstates,
            neighbor_partial,
            log_scale);

        Node *true_child = (parent_of[neighbor->id] == node->id) ? neighbor : node;
        const double *transition =
            &transition_cache[(size_t)true_child->id * nstates * nstates];
        if (true_child == neighbor) {
            // Downward edge: contribution(x) = sum_y P(x, y) message(y).
            for (int node_state = 0; node_state < nstates; node_state++) {
                double contribution = 0.0;
                for (int child_state = 0; child_state < nstates; child_state++) {
                    contribution += transition[node_state * nstates + child_state] *
                        neighbor_partial[child_state];
                }
                partial[node_state] *= contribution;
            }
        } else {
            // Upward edge: node is the true child, so the message from above
            // transports along the transposed direction:
            // contribution(x) = sum_y message(y) P(y, x).
            for (int node_state = 0; node_state < nstates; node_state++) {
                double contribution = 0.0;
                for (int parent_state = 0; parent_state < nstates; parent_state++) {
                    contribution += neighbor_partial[parent_state] *
                        transition[parent_state * nstates + node_state];
                }
                partial[node_state] *= contribution;
            }
        }
    }

    // Same rescaling as the reversible recursion; any positive scale cancels
    // in the normalised posterior downstream.
    double partial_sum = 0.0;
    for (int state = 0; state < nstates; state++)
        partial_sum += partial[state];
    if (partial_sum > 0.0 && isfinite(partial_sum)) {
        if (log_scale != nullptr)
            *log_scale += log(partial_sum);
        for (int state = 0; state < nstates; state++)
            partial[state] /= partial_sum;
    }
}

/** Position of `target` in `node`'s neighbour list; -1 when not adjacent. */
static int neighborIndex(Node *node, Node *target) {
    for (size_t i = 0; i < node->neighbors.size(); i++)
        if (node->neighbors[i]->node == target)
            return (int)i;
    return -1;
}

/** Same rescaling the recursion applies; any positive scale cancels later. */
static inline void rescalePartial(double *partial, int nstates, double *log_scale = nullptr) {
    double partial_sum = 0.0;
    for (int state = 0; state < nstates; state++)
        partial_sum += partial[state];
    if (partial_sum > 0.0 && isfinite(partial_sum)) {
        if (log_scale != nullptr)
            *log_scale += log(partial_sum);
        for (int state = 0; state < nstates; state++)
            partial[state] /= partial_sum;
    }
}

size_t directionalCacheBytes(
    size_t node_bound, int max_degree, size_t nptn, int nstates, bool track_scales)
{
    return node_bound * (size_t)max_degree * nptn * (nstates + (track_scales ? 1 : 0)) * sizeof(double);
}

/**
 * Post-order half of the two-pass sweep: fills msg(node -> dad) for `node` and
 * everything below it. Identical arithmetic to computeOrdinaryPartial(), just
 * recorded instead of discarded.
 */
static void upwardMessages(
    Node *node,
    Node *dad,
    size_t ptn,
    Alignment *aln,
    ModelSubst *model,
    int nstates,
    const vector<double> &transition_cache,
    SatuTeDirectionalCache &cache)
{
    double *out = cache.message(node->id, neighborIndex(node, dad), ptn);

    if (node->isLeaf()) {
        StateType state = (node->id >= 0 && node->id < (int)aln->getNSeq())
            ? aln->at(ptn)[node->id]
            : aln->STATE_UNKNOWN;
        model->computeTipLikelihood(state, out);
        return;
    }

    FOR_NEIGHBOR_IT(node, dad, it)
        upwardMessages((*it)->node, node, ptn, aln, model, nstates, transition_cache, cache);

    for (int state = 0; state < nstates; state++)
        out[state] = 1.0;
    double scale = 0.0;
    FOR_NEIGHBOR_IT(node, dad, it) {
        Node *child = (*it)->node;
        const double *child_msg =
            cache.message(child->id, neighborIndex(child, node), ptn);
        if (!cache.log_scales.empty())
            scale += cache.logScale(child->id, neighborIndex(child, node), ptn);
        // Transition matrices are keyed by the node the edge descends into, so
        // under this rooting the key for edge (node, child) is child's id.
        const double *transition =
            &transition_cache[(size_t)child->id * nstates * nstates];
        for (int parent_state = 0; parent_state < nstates; parent_state++) {
            double contribution = 0.0;
            for (int child_state = 0; child_state < nstates; child_state++) {
                contribution += transition[parent_state * nstates + child_state] *
                    child_msg[child_state];
            }
            out[parent_state] *= contribution;
        }
    }
    rescalePartial(out, nstates, cache.log_scales.empty() ? nullptr : &scale);
    if (!cache.log_scales.empty())
        cache.logScale(node->id, neighborIndex(node, dad), ptn) = scale;
}

/**
 * Pre-order half: with every msg(neighbour -> node) already known, emit
 * msg(node -> child) for each child. The message toward `dad` is skipped
 * because the post-order pass already wrote it.
 *
 * Each outgoing message is the product of all incoming messages except the one
 * it answers, so prefix and suffix products give all of them in one sweep
 * rather than one product per outgoing edge.
 */
static void downwardMessages(
    Node *node,
    Node *dad,
    size_t ptn,
    int nstates,
    const vector<double> &transition_cache,
    SatuTeDirectionalCache &cache,
    vector<double> &incoming,
    vector<double> &prefix,
    vector<double> &suffix)
{
    if (node->isLeaf())
        return;

    const int degree = (int)node->neighbors.size();
    double total_scale = 0.0;
    for (int i = 0; i < degree; i++) {
        Node *other = node->neighbors[i]->node;
        // Under the global rooting the deeper endpoint owns the matrix: that is
        // `node` itself for the edge back to dad, and the neighbour otherwise.
        const int key = (other == dad) ? node->id : other->id;
        const double *in_msg = cache.message(other->id, neighborIndex(other, node), ptn);
        if (!cache.log_scales.empty())
            total_scale += cache.logScale(other->id, neighborIndex(other, node), ptn);
        const double *transition = &transition_cache[(size_t)key * nstates * nstates];
        double *dst = &incoming[(size_t)i * nstates];
        for (int parent_state = 0; parent_state < nstates; parent_state++) {
            double contribution = 0.0;
            for (int child_state = 0; child_state < nstates; child_state++) {
                contribution += transition[parent_state * nstates + child_state] *
                    in_msg[child_state];
            }
            dst[parent_state] = contribution;
        }
    }

    for (int state = 0; state < nstates; state++) {
        prefix[state] = 1.0;
        suffix[(size_t)degree * nstates + state] = 1.0;
    }
    for (int i = 0; i < degree; i++) {
        for (int state = 0; state < nstates; state++) {
            prefix[(size_t)(i + 1) * nstates + state] =
                prefix[(size_t)i * nstates + state] * incoming[(size_t)i * nstates + state];
        }
    }
    for (int i = degree - 1; i >= 0; i--) {
        for (int state = 0; state < nstates; state++) {
            suffix[(size_t)i * nstates + state] =
                suffix[(size_t)(i + 1) * nstates + state] * incoming[(size_t)i * nstates + state];
        }
    }

    for (int j = 0; j < degree; j++) {
        if (node->neighbors[j]->node == dad)
            continue;
        double *out = cache.message(node->id, j, ptn);
        for (int state = 0; state < nstates; state++) {
            out[state] = prefix[(size_t)j * nstates + state] *
                suffix[(size_t)(j + 1) * nstates + state];
        }
        double scale = 0.0;
        if (!cache.log_scales.empty()) {
            Node *other = node->neighbors[j]->node;
            scale = total_scale - cache.logScale(other->id, neighborIndex(other, node), ptn);
        }
        rescalePartial(out, nstates, cache.log_scales.empty() ? nullptr : &scale);
        if (!cache.log_scales.empty())
            cache.logScale(node->id, j, ptn) = scale;
    }

    // Safe to clobber the scratch below: every outgoing message is already
    // written into the cache.
    FOR_NEIGHBOR_IT(node, dad, it) {
        downwardMessages((*it)->node, node, ptn, nstates, transition_cache, cache,
                         incoming, prefix, suffix);
    }
}

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
    bool track_scales)
{
    cache.ready = false;
    if (root == nullptr || root->isLeaf() || max_degree <= 0)
        return false;

    const size_t nptn = aln->size();
    cache.max_degree = max_degree;
    cache.nstates = nstates;
    cache.nptn = nptn;
    cache.dir_stride = nptn * (size_t)nstates;
    cache.data.assign(node_bound * (size_t)max_degree * cache.dir_stride, 0.0);
    cache.log_scales.clear();
    if (track_scales)
        cache.log_scales.assign(node_bound * (size_t)max_degree * nptn, 0.0);

    // One matrix per edge, filled from the same rooting the sweeps use.
    vector<double> transition_cache(node_bound * (size_t)(nstates * nstates), 0.0);
    FOR_NEIGHBOR_IT(root, nullptr, it) {
        ensureTransitionMatrix(*it, model, nstates, rate_multiplier, transition_cache);
        prepareTransitionMatrices((*it)->node, root, model, nstates, rate_multiplier,
                                  transition_cache);
    }

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = (num_threads > 0) ? num_threads : omp_get_max_threads();
#endif

    const int nptn_int = (int)nptn;
#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads)
#endif
    {
        vector<double> incoming((size_t)max_degree * nstates, 0.0);
        vector<double> prefix((size_t)(max_degree + 1) * nstates, 0.0);
        vector<double> suffix((size_t)(max_degree + 1) * nstates, 0.0);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (int ptn = 0; ptn < nptn_int; ptn++) {
            FOR_NEIGHBOR_IT(root, nullptr, it) {
                upwardMessages((*it)->node, root, (size_t)ptn, aln, model, nstates,
                               transition_cache, cache);
            }
            downwardMessages(root, nullptr, (size_t)ptn, nstates, transition_cache,
                             cache, incoming, prefix, suffix);
        }
    }

    cache.ready = true;
    return true;
}

static void finalizeSatuTeResult(SatuTeBranchResult &result) {
    if (result.valid_sites <= 0) {
        result.decision = "no_sites";
        result.decision_bonferroni = "no_sites";
        return;
    }
    if (!(result.variance > 0.0) || !isfinite(result.coherence) || !isfinite(result.variance)) {
        result.decision = "undefined_variance";
        result.decision_bonferroni = "undefined_variance";
        return;
    }

    result.se = sqrt(result.variance / result.valid_sites);
    if (result.se > 0.0 && isfinite(result.se)) {
        result.z_score = result.coherence / result.se;
        result.p_value = gsl_cdf_ugaussian_Q(result.z_score);
        result.decision = (result.p_value <= result.alpha) ? "informative" : "saturated";
        result.decision_bonferroni = (result.p_value <= result.alpha_adjusted) ? "informative" : "saturated";
    }
}

void jacobiEigenSymmetric(vector<double> &a, int n, vector<double> &eigenvectors)
{
    eigenvectors.assign((size_t)n * n, 0.0);
    for (int i = 0; i < n; i++)
        eigenvectors[(size_t)i * n + i] = 1.0;

    const int max_sweeps = 100;
    for (int sweep = 0; sweep < max_sweeps; sweep++) {
        double off_diagonal = 0.0;
        for (int p = 0; p < n; p++)
            for (int q = p + 1; q < n; q++)
                off_diagonal += a[(size_t)p * n + q] * a[(size_t)p * n + q];
        if (off_diagonal < 1e-30)
            break;

        for (int p = 0; p < n; p++) {
            for (int q = p + 1; q < n; q++) {
                if (fabs(a[(size_t)p * n + q]) < 1e-300)
                    continue;
                double theta = (a[(size_t)q * n + q] - a[(size_t)p * n + p]) /
                    (2.0 * a[(size_t)p * n + q]);
                double t = ((theta >= 0.0) ? 1.0 : -1.0) / (fabs(theta) + sqrt(theta * theta + 1.0));
                double c = 1.0 / sqrt(t * t + 1.0);
                double s = t * c;
                double app = a[(size_t)p * n + p];
                double aqq = a[(size_t)q * n + q];
                double apq = a[(size_t)p * n + q];

                a[(size_t)p * n + p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
                a[(size_t)q * n + q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
                a[(size_t)p * n + q] = 0.0;
                a[(size_t)q * n + p] = 0.0;

                for (int i = 0; i < n; i++) {
                    if (i == p || i == q)
                        continue;
                    double aip = a[(size_t)i * n + p];
                    double aiq = a[(size_t)i * n + q];
                    a[(size_t)i * n + p] = c * aip - s * aiq;
                    a[(size_t)p * n + i] = a[(size_t)i * n + p];
                    a[(size_t)i * n + q] = s * aip + c * aiq;
                    a[(size_t)q * n + i] = a[(size_t)i * n + q];
                }
                for (int i = 0; i < n; i++) {
                    double vip = eigenvectors[(size_t)i * n + p];
                    double viq = eigenvectors[(size_t)i * n + q];
                    eigenvectors[(size_t)i * n + p] = c * vip - s * viq;
                    eigenvectors[(size_t)i * n + q] = s * vip + c * viq;
                }
            }
        }
    }
}

/**
 * The model mean is M a, so M+ times that mean is the projection M+ M a.
 * Retain persistence weights on identifiable directions; do not invert a
 * second time by mistaking a itself for the model mean.
 */
static vector<double> solveEnhancedWeights(const vector<double> &noise_matrix, const vector<double> &signal, size_t n)
{
    if (n == 1) {
        double m = noise_matrix[0];
        double tolerance = fabs(m) * 1e-8;
        return vector<double>(1, (m > tolerance) ? signal[0] : 0.0);
    }

    vector<double> a = noise_matrix;
    vector<double> eigenvectors;
    jacobiEigenSymmetric(a, (int)n, eigenvectors);

    vector<double> eigenvalues(n);
    double max_eigenvalue = 0.0;
    for (size_t i = 0; i < n; i++) {
        eigenvalues[i] = a[i * n + i];
        max_eigenvalue = max(max_eigenvalue, eigenvalues[i]);
    }
    double tolerance = max_eigenvalue * 1e-8;

    vector<double> weights(n, 0.0);
    for (size_t k = 0; k < n; k++) {
        if (eigenvalues[k] <= tolerance)
            continue;
        double projection = 0.0;
        for (size_t i = 0; i < n; i++)
            projection += eigenvectors[i * n + k] * signal[i];
        double scale = projection;
        for (size_t i = 0; i < n; i++)
            weights[i] += scale * eigenvectors[i * n + k];
    }
    return weights;
}

static void addMixtureCategory(
    SatuTeMixtureAccumulator &mixture, size_t ptn,
    double log_density, double category_score)
{
    double previous = mixture.log_null[ptn];
    double largest = max(previous, log_density);
    double combined = largest + log(exp(previous - largest) + exp(log_density - largest));
    mixture.score[ptn] = mixture.score[ptn] * exp(previous - combined)
        + category_score * exp(log_density - combined);
    mixture.log_null[ptn] = combined;
}

static SatuTeBranchResult computeCategoryStatisticImpl(
    Node *left_node,
    Node *left_dad,
    Node *right_node,
    Node *right_dad,
    Alignment *aln,
    ModelSubst *model,
    int nstates,
    const double *eval,
    const SatuTeFormulaSpec &formula,
    const vector<double> &mode_weights,
    const IntVector *pattern_cat,
    int required_pattern_cat,
    double rate_multiplier,
    int num_threads,
    bool enhanced_weights,
    const SatuTeBranchResult &base,
    const SatuTeDirectionalCache *directional_cache = nullptr,
    vector<double> *out_pattern_coherence = nullptr,
    const SatuTeDirectedContext *directed = nullptr,
    SatuTeMixtureAccumulator *mixture = nullptr)
{
    SatuTeBranchResult result = base;
    size_t nmodes = formula.modes.size();
    vector<double> left_second(nmodes * nmodes, 0.0);
    vector<double> right_second(nmodes * nmodes, 0.0);
    // First moments, needed to centre the statistic. The raw cross-moment
    // sum_k w_k L_k R_k is only a valid independence test when E[L_k]=E[R_k]=0,
    // which holds when the fitted model's eigenvectors are orthogonal to its
    // own stationary distribution. Under a misspecified model both means drift
    // off zero and E[C] = sum_k w_k E[L_k] E[R_k] != 0, a bias that grows as
    // sqrt(n) in the z-score. Centring turns the statistic into a covariance,
    // which is zero under independence whatever the marginal means are.
    vector<double> left_first(nmodes, 0.0);
    vector<double> right_first(nmodes, 0.0);
    // Per-mode, unweighted: sum_site freq * left_factor[i] * right_factor[i].
    // Weighting is linear, so summing per mode and combining with a weight
    // vector afterward equals summing pre-weighted per site, but deferring it
    // is what lets eigenvalue_weighted solve for its weights from the same
    // covariance matrices this loop already produces, without a second pass
    // over the alignment.
    vector<double> coherence_vec(nmodes, 0.0);
    vector<double> state_freq(nstates, 0.0);
    // One matrix slot per node id; 2*nseq is a safe upper bound on node count.
    vector<double> transition_cache(
        (2 * aln->getNSeq() + 2) * (size_t)(nstates * nstates), 0.0);
    double *eigenvectors = model->getEigenvectors();

    model->getStateFrequency(&state_freq[0]);

    // Directed non-reversible evaluation recomputes both side vectors with the
    // direction-aware recursion; the shared reversible message cache does not
    // apply to it.
    const bool use_directed = (directed != nullptr);

    // When the caller supplied the shared directional messages for this rate
    // category, both partials this branch needs are already computed and the
    // per-branch subtree recursion (and its transition matrices) is skipped
    // entirely. That is what removes the O(n^2)-in-taxa cost.
    const bool use_cache =
        !use_directed && (directional_cache != nullptr) && directional_cache->ready;
    const int left_message_index = use_cache ? neighborIndex(left_node, left_dad) : -1;
    const int right_message_index = use_cache ? neighborIndex(right_node, right_dad) : -1;

    vector<char> transition_filled;
    if (use_directed) {
        transition_filled.assign(2 * aln->getNSeq() + 2, 0);
        prepareDirectedTransitionMatrices(
            left_node, left_dad, model, nstates, rate_multiplier,
            *directed->parent_of, transition_cache, transition_filled);
        prepareDirectedTransitionMatrices(
            right_node, right_dad, model, nstates, rate_multiplier,
            *directed->parent_of, transition_cache, transition_filled);
    } else if (!use_cache) {
        // Hoisted out of the pattern loop: after this the cache is read-only,
        // which is what makes the loop below safe to share between threads.
        prepareTransitionMatrices(
            left_node, left_dad, model, nstates, rate_multiplier, transition_cache);
        prepareTransitionMatrices(
            right_node, right_dad, model, nstates, rate_multiplier, transition_cache);
    }

    const int nptn = (int)aln->size();

    // Per-pattern factors, retained only when a caller asked for the per-site
    // export. The centred, weighted per-pattern contribution cannot be formed
    // inside the loop because both the marginal means and (for
    // eigenvalue_weighted with enhanced weights) the weights themselves are
    // only known after the whole alignment has been summed -- so the raw
    // factors are kept and combined in a cheap second pass below.
    // schedule(static) gives each pattern to exactly one thread, so these
    // shared arrays are written race-free without extra synchronisation.
    const bool want_pattern_coherence = (out_pattern_coherence != nullptr);
    vector<double> pattern_left_factor;
    vector<double> pattern_right_factor;
    vector<char> pattern_counted;
    if (want_pattern_coherence) {
        pattern_left_factor.assign((size_t)nptn * nmodes, 0.0);
        pattern_right_factor.assign((size_t)nptn * nmodes, 0.0);
        pattern_counted.assign((size_t)nptn, 0);
        out_pattern_coherence->assign((size_t)nptn, 0.0);
    }

    // The mixture score is centred for the same reason the homogeneous one is:
    // the raw cross-moment only tests independence when E[L_k]=E[R_k]=0, which
    // a misspecified model breaks. The marginal means are over every valid site
    // rather than this category's assigned subset, so like the export above the
    // factors are kept and the score formed in a second pass.
    const bool want_mixture = (mixture != nullptr && mixture->proportion > 0.0);
    vector<double> mixture_left_factor;
    vector<double> mixture_right_factor;
    vector<double> mixture_log_density;
    vector<char> mixture_counted;
    if (want_mixture) {
        mixture_left_factor.assign((size_t)nptn * nmodes, 0.0);
        mixture_right_factor.assign((size_t)nptn * nmodes, 0.0);
        mixture_log_density.assign((size_t)nptn, 0.0);
        mixture_counted.assign((size_t)nptn, 0);
    }

    // Patterns are independent, so the loop parallelises over them. Each thread
    // accumulates into its own slot and the slots are combined in thread order
    // after the region; together with schedule(static) that makes results bit
    // reproducible for a given thread count, rather than depending on which
    // thread happens to finish first. Per-thread slots also avoid OpenMP 4.5
    // array reductions, which MSVC does not support.
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = (num_threads > 0) ? num_threads : omp_get_max_threads();
#endif
    const size_t second_size = nmodes * nmodes;
    vector<double> thread_valid(nthreads, 0.0);
    vector<int> thread_rate_category_sites(nthreads, 0);
    vector<int> thread_skipped(nthreads, 0);
    vector<double> thread_left_second(nthreads * second_size, 0.0);
    vector<double> thread_right_second(nthreads * second_size, 0.0);
    vector<double> thread_left_first(nthreads * nmodes, 0.0);
    vector<double> thread_right_first(nthreads * nmodes, 0.0);
    vector<double> thread_coherence_vec(nthreads * nmodes, 0.0);
    vector<double> thread_mixture_left_first(nthreads * nmodes, 0.0);
    vector<double> thread_mixture_right_first(nthreads * nmodes, 0.0);
    vector<double> thread_mixture_sites(nthreads, 0.0);

#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads)
#endif
    {
    int thread_id = 0;
#ifdef _OPENMP
    thread_id = omp_get_thread_num();
#endif
    vector<double> left_factor(nmodes, 0.0);
    vector<double> right_factor(nmodes, 0.0);
    vector<double> left_likelihood(nstates, 0.0);
    vector<double> right_likelihood(nstates, 0.0);
    vector<double> left_posterior(nstates, 0.0);
    vector<double> right_posterior(nstates, 0.0);
    // Depth of the recursion is bounded by the number of leaves; one nstates
    // block per level, allocated once per thread rather than per node visit.
    vector<double> scratch((aln->getNSeq() + 2) * (size_t)nstates, 0.0);
    vector<double> local_left_second(second_size, 0.0);
    vector<double> local_right_second(second_size, 0.0);
    vector<double> local_left_first(nmodes, 0.0);
    vector<double> local_right_first(nmodes, 0.0);
    vector<double> local_coherence_vec(nmodes, 0.0);
    vector<double> local_mixture_left_first(nmodes, 0.0);
    vector<double> local_mixture_right_first(nmodes, 0.0);
    double local_mixture_sites = 0.0;
    double local_valid_sites = 0.0;
    int local_rate_category_sites = 0;
    int local_skipped_sites = 0;

    // Static: rate categories are scattered across patterns rather than
    // clustered, so fixed chunks stay balanced, and the mapping from pattern to
    // thread is reproducible.
#ifdef _OPENMP
#pragma omp for schedule(static) nowait
#endif
    for (int ptn = 0; ptn < nptn; ptn++) {
        bool assigned = true;
        if (pattern_cat != nullptr && required_pattern_cat >= 0) {
            assigned = (size_t)ptn < pattern_cat->size() && (*pattern_cat)[ptn] == required_pattern_cat;
        }
        if (!assigned && !want_mixture)
            continue;

        int freq = aln->at(ptn).frequency;
        if (freq <= 0)
            continue;
        if (assigned)
            local_rate_category_sites += freq;

        const double *left_lh;
        const double *right_lh;
        double left_log_scale = 0.0;
        double right_log_scale = 0.0;
        if (use_cache) {
            left_lh = directional_cache->message(left_node->id, left_message_index, ptn);
            right_lh = directional_cache->message(right_node->id, right_message_index, ptn);
            if (want_mixture) {
                left_log_scale = directional_cache->logScale(left_node->id, left_message_index, ptn);
                right_log_scale = directional_cache->logScale(right_node->id, right_message_index, ptn);
            }
        } else if (use_directed) {
            computeDirectedPartial(
                left_node,
                left_dad,
                ptn,
                aln,
                model,
                nstates,
                *directed->parent_of,
                directed->root_node_id,
                &state_freq[0],
                transition_cache,
                &scratch[0],
                &left_likelihood[0], !want_mixture ? nullptr : &left_log_scale);
            computeDirectedPartial(
                right_node,
                right_dad,
                ptn,
                aln,
                model,
                nstates,
                *directed->parent_of,
                directed->root_node_id,
                &state_freq[0],
                transition_cache,
                &scratch[0],
                &right_likelihood[0], !want_mixture ? nullptr : &right_log_scale);
            left_lh = &left_likelihood[0];
            right_lh = &right_likelihood[0];
        } else {
            computeOrdinaryPartial(
                left_node,
                left_dad,
                ptn,
                aln,
                model,
                nstates,
                transition_cache,
                &scratch[0],
                &left_likelihood[0], !want_mixture ? nullptr : &left_log_scale);
            computeOrdinaryPartial(
                right_node,
                right_dad,
                ptn,
                aln,
                model,
                nstates,
                transition_cache,
                &scratch[0],
                &right_likelihood[0], !want_mixture ? nullptr : &right_log_scale);
            left_lh = &left_likelihood[0];
            right_lh = &right_likelihood[0];
        }

        double left_sum = 0.0;
        double right_sum = 0.0;
        for (int state = 0; state < nstates; state++) {
            // Directed rootward vectors already carry the stationary prior
            // from the root leaf, so only the tipward side multiplies pi.
            left_posterior[state] = use_directed
                ? left_lh[state]
                : left_lh[state] * state_freq[state];
            right_posterior[state] = right_lh[state] * state_freq[state];
            left_sum += left_posterior[state];
            right_sum += right_posterior[state];
        }

        if (!isfinite(left_sum) || !isfinite(right_sum) ||
            fabs(left_sum) <= DBL_MIN || fabs(right_sum) <= DBL_MIN) {
            if (assigned)
                local_skipped_sites += freq;
            if (want_mixture && (!isfinite(left_sum) || !isfinite(right_sum)))
                mixture->invalid[ptn] = 1;
            continue;
        }
        for (int state = 0; state < nstates; state++) {
            left_posterior[state] /= left_sum;
            right_posterior[state] /= right_sum;
        }

        bool valid = true;
        for (size_t i = 0; i < nmodes; i++) {
            left_factor[i] = 0.0;
            right_factor[i] = 0.0;
            if (use_directed) {
                // Non-reversible: each side has its own singular basis.
                for (int state = 0; state < nstates; state++) {
                    left_factor[i] += directed->left_projection[i * nstates + state] *
                        left_posterior[state];
                    right_factor[i] += directed->right_projection[i * nstates + state] *
                        right_posterior[state];
                }
            } else {
                for (int state = 0; state < nstates; state++) {
                    double eigenvector_value = eigenvectors[state * nstates + formula.modes[i]];
                    left_factor[i] += eigenvector_value * left_posterior[state];
                    right_factor[i] += eigenvector_value * right_posterior[state];
                }
            }
            if (!isfinite(left_factor[i]) || !isfinite(right_factor[i])) {
                valid = false;
                break;
            }
        }
        if (!valid) {
            if (assigned)
                local_skipped_sites += freq;
            if (want_mixture)
                mixture->invalid[ptn] = 1;
            continue;
        }

        if (want_mixture) {
            // A zero-rate category keeps its joint probability in the null but
            // carries no branch signal, so it stores no factors and scores 0.
            const bool has_signal = (rate_multiplier != 0.0);
            double null_ratio = 1.0;
            if (!has_signal) {
                null_ratio = 0.0;
                for (int state = 0; state < nstates; state++)
                    if (state_freq[state] > 0.0)
                        null_ratio += left_posterior[state] * right_posterior[state] / state_freq[state];
            }
            if (null_ratio > 0.0) {
                double log_density = log(mixture->proportion) + log(left_sum) + log(right_sum)
                    + left_log_scale + right_log_scale + log(null_ratio);
                if (isfinite(log_density)) {
                    mixture_log_density[ptn] = log_density;
                    mixture_counted[ptn] = 1;
                    if (has_signal) {
                        for (size_t i = 0; i < nmodes; i++) {
                            mixture_left_factor[(size_t)ptn * nmodes + i] = left_factor[i];
                            mixture_right_factor[(size_t)ptn * nmodes + i] = right_factor[i];
                            local_mixture_left_first[i] += freq * left_factor[i];
                            local_mixture_right_first[i] += freq * right_factor[i];
                        }
                        local_mixture_sites += freq;
                    }
                } else {
                    mixture->invalid[ptn] = 1;
                }
            }
        }
        if (!assigned)
            continue;

        local_valid_sites += freq;
        if (want_pattern_coherence) {
            pattern_counted[ptn] = 1;
            for (size_t i = 0; i < nmodes; i++) {
                pattern_left_factor[(size_t)ptn * nmodes + i] = left_factor[i];
                pattern_right_factor[(size_t)ptn * nmodes + i] = right_factor[i];
            }
        }
        for (size_t i = 0; i < nmodes; i++) {
            local_left_first[i] += freq * left_factor[i];
            local_right_first[i] += freq * right_factor[i];
            local_coherence_vec[i] += freq * left_factor[i] * right_factor[i];
        }
        for (size_t i = 0; i < nmodes; i++) {
            for (size_t j = 0; j < nmodes; j++) {
                size_t idx = i * nmodes + j;
                local_left_second[idx] += freq * left_factor[i] * left_factor[j];
                local_right_second[idx] += freq * right_factor[i] * right_factor[j];
            }
        }
    }

    thread_valid[thread_id] = local_valid_sites;
    thread_rate_category_sites[thread_id] = local_rate_category_sites;
    thread_skipped[thread_id] = local_skipped_sites;
    for (size_t idx = 0; idx < second_size; idx++) {
        thread_left_second[thread_id * second_size + idx] = local_left_second[idx];
        thread_right_second[thread_id * second_size + idx] = local_right_second[idx];
    }
    for (size_t i = 0; i < nmodes; i++) {
        thread_left_first[thread_id * nmodes + i] = local_left_first[i];
        thread_right_first[thread_id * nmodes + i] = local_right_first[i];
        thread_coherence_vec[thread_id * nmodes + i] = local_coherence_vec[i];
        thread_mixture_left_first[thread_id * nmodes + i] = local_mixture_left_first[i];
        thread_mixture_right_first[thread_id * nmodes + i] = local_mixture_right_first[i];
    }
    thread_mixture_sites[thread_id] = local_mixture_sites;
    } // omp parallel

    // Fixed combination order, so the totals do not depend on thread timing.
    double valid_sites = 0.0;
    int rate_category_sites = 0;
    int skipped_sites = 0;
    for (int thread = 0; thread < nthreads; thread++) {
        valid_sites += thread_valid[thread];
        rate_category_sites += thread_rate_category_sites[thread];
        skipped_sites += thread_skipped[thread];
        for (size_t idx = 0; idx < second_size; idx++) {
            left_second[idx] += thread_left_second[thread * second_size + idx];
            right_second[idx] += thread_right_second[thread * second_size + idx];
        }
        for (size_t i = 0; i < nmodes; i++) {
            left_first[i] += thread_left_first[thread * nmodes + i];
            right_first[i] += thread_right_first[thread * nmodes + i];
            coherence_vec[i] += thread_coherence_vec[thread * nmodes + i];
        }
    }

    // Before the early return below: a category can own no assigned pattern
    // yet still owe the mixture its share of every site.
    if (want_mixture) {
        double mixture_sites = 0.0;
        vector<double> mixture_left_mean(nmodes, 0.0);
        vector<double> mixture_right_mean(nmodes, 0.0);
        for (int thread = 0; thread < nthreads; thread++) {
            mixture_sites += thread_mixture_sites[thread];
            for (size_t i = 0; i < nmodes; i++) {
                mixture_left_mean[i] += thread_mixture_left_first[thread * nmodes + i];
                mixture_right_mean[i] += thread_mixture_right_first[thread * nmodes + i];
            }
        }
        if (mixture_sites > 0.0) {
            for (size_t i = 0; i < nmodes; i++) {
                mixture_left_mean[i] /= mixture_sites;
                mixture_right_mean[i] /= mixture_sites;
            }
        }
        // Serial and in pattern order: the streaming combination is not
        // associative, so this keeps results reproducible across thread counts.
        for (int ptn = 0; ptn < nptn; ptn++) {
            if (!mixture_counted[ptn])
                continue;
            double category_score = 0.0;
            for (size_t i = 0; i < nmodes; i++) {
                double centred_left =
                    mixture_left_factor[(size_t)ptn * nmodes + i] - mixture_left_mean[i];
                double centred_right =
                    mixture_right_factor[(size_t)ptn * nmodes + i] - mixture_right_mean[i];
                category_score += mode_weights[i] * centred_left * centred_right;
            }
            category_score *= mixture->weight_scale;
            if (isfinite(category_score))
                addMixtureCategory(*mixture, (size_t)ptn, mixture_log_density[ptn], category_score);
            else
                mixture->invalid[ptn] = 1;
        }
    }

    result.valid_sites = (int)valid_sites;
    result.rate_category_sites = rate_category_sites;
    result.skipped_sites = skipped_sites;
    if (valid_sites <= 0.0) {
        finalizeSatuTeResult(result);
        return result;
    }

    // Centre with the observed marginal means. Under a correctly specified
    // model these are ~0 and this is a no-op; under misspecification they carry
    // the bias that would otherwise be read as branch signal.
    vector<double> left_mean(nmodes, 0.0);
    vector<double> right_mean(nmodes, 0.0);
    for (size_t i = 0; i < nmodes; i++) {
        left_mean[i] = left_first[i] / valid_sites;
        right_mean[i] = right_first[i] / valid_sites;
    }

    // Centred second moments -- Cov_L(i,j), Cov_R(i,j). Needed for the
    // variance regardless of formula, and, for eigenvalue_weighted only, to
    // build the noise-covariance matrix the enhanced weights solve against.
    vector<double> left_cov(nmodes * nmodes, 0.0);
    vector<double> right_cov(nmodes * nmodes, 0.0);
    for (size_t i = 0; i < nmodes; i++) {
        for (size_t j = 0; j < nmodes; j++) {
            size_t idx = i * nmodes + j;
            left_cov[idx] = left_second[idx] / valid_sites - left_mean[i] * left_mean[j];
            right_cov[idx] = right_second[idx] / valid_sites - right_mean[i] * right_mean[j];
        }
    }

    // The optional projection removes unidentifiable directions of M a.
    vector<double> active_weights = mode_weights;
    if (formula.eigenvalue_weighted && enhanced_weights) {
        vector<double> noise_matrix(nmodes * nmodes, 0.0);
        for (size_t idx = 0; idx < nmodes * nmodes; idx++)
            noise_matrix[idx] = left_cov[idx] * right_cov[idx];
        active_weights = solveEnhancedWeights(noise_matrix, mode_weights, nmodes);
    }

    double centred_coherence = 0.0;
    for (size_t i = 0; i < nmodes; i++)
        centred_coherence += active_weights[i] * (coherence_vec[i] / valid_sites - left_mean[i] * right_mean[i]);
    result.coherence = centred_coherence;

    // Second pass for the per-site export: now that the means and the actually
    // applied weights (enhanced or not) are known, each retained pattern's
    // centred weighted contribution can be formed. Frequency weighted and
    // divided by valid_sites these sum to centred_coherence above.
    if (want_pattern_coherence) {
        for (int ptn = 0; ptn < nptn; ptn++) {
            if (!pattern_counted[ptn])
                continue;
            double contribution = 0.0;
            for (size_t i = 0; i < nmodes; i++) {
                double centred_left = pattern_left_factor[(size_t)ptn * nmodes + i] - left_mean[i];
                double centred_right = pattern_right_factor[(size_t)ptn * nmodes + i] - right_mean[i];
                contribution += active_weights[i] * centred_left * centred_right;
            }
            (*out_pattern_coherence)[ptn] = contribution;
        }
    }

    // Matching centred second moments, so the variance describes the same
    // centred quantity the numerator now measures.
    result.variance = 0.0;
    for (size_t i = 0; i < nmodes; i++) {
        for (size_t j = 0; j < nmodes; j++) {
            size_t idx = i * nmodes + j;
            result.variance += active_weights[i] * active_weights[j] * left_cov[idx] * right_cov[idx];
        }
    }
    // Report what was actually applied. For dominant this is identical to
    // the pre-solve weights already formatted by the caller; for
    // eigenvalue_weighted it is now data-dependent and only known here.
    result.weights = format::formatVector(active_weights);
    finalizeSatuTeResult(result);
    if (mixture != nullptr) {
        result.p_value = numeric_limits<double>::quiet_NaN();
        result.decision = "diagnostic";
        result.decision_bonferroni = "not_tested";
    }
    return result;
}

SatuTeBranchResult computeMixtureStatistic(
    const SatuTeMixtureAccumulator &mixture, Alignment *aln,
    const SatuTeBranchResult &base, vector<double> *out_pattern_coherence)
{
    SatuTeBranchResult result = base;
    result.coherence = 0.0;
    if (out_pattern_coherence != nullptr)
        out_pattern_coherence->assign(aln->size(), 0.0);
    for (size_t ptn = 0; ptn < aln->size(); ptn++) {
        int freq = aln->at(ptn).frequency;
        result.rate_category_sites += freq;
        if (mixture.invalid[ptn] || !isfinite(mixture.log_null[ptn]) || !isfinite(mixture.score[ptn])) {
            result.skipped_sites += freq;
            continue;
        }
        result.valid_sites += freq;
        result.coherence += freq * mixture.score[ptn];
        if (out_pattern_coherence != nullptr)
            (*out_pattern_coherence)[ptn] = mixture.score[ptn];
    }
    if (result.valid_sites > 0)
        result.coherence /= result.valid_sites;
    result.variance = 0.0;
    for (size_t ptn = 0; ptn < aln->size(); ptn++) {
        if (mixture.invalid[ptn] || !isfinite(mixture.log_null[ptn]) || !isfinite(mixture.score[ptn]))
            continue;
        double difference = mixture.score[ptn] - result.coherence;
        result.variance += aln->at(ptn).frequency * difference * difference;
    }
    if (result.valid_sites > 1)
        result.variance /= result.valid_sites - 1;
    finalizeSatuTeResult(result);
    return result;
}

void applySeparateFormulaFdr(
    vector<SatuTeBranchResult> &results,
    double fdr_level)
{
    const char *formula_families[] = {"dominant", "eigenvalue_weighted"};

    for (size_t result_index = 0; result_index < results.size(); result_index++) {
        results[result_index].fdr_by = numeric_limits<double>::quiet_NaN();
        results[result_index].decision_fdr =
            (results[result_index].rate_category == "pooled")
                ? results[result_index].decision
                : "not_tested";
    }

    for (size_t family_index = 0; family_index < 2; family_index++) {
        vector<pair<double, size_t> > ordered;
        for (size_t result_index = 0; result_index < results.size(); result_index++) {
            const SatuTeBranchResult &result = results[result_index];
            if (result.formula != formula_families[family_index] ||
                result.rate_category != "pooled")
                continue;
            // Undefined tests remain in the planned branch family as p=1 so
            // they cannot be rejected or reduce the multiplicity penalty.
            double family_p_value =
                (isfinite(result.p_value) && result.p_value >= 0.0 && result.p_value <= 1.0)
                    ? result.p_value
                    : 1.0;
            ordered.push_back(make_pair(family_p_value, result_index));
        }

        stable_sort(
            ordered.begin(),
            ordered.end(),
            [](const pair<double, size_t> &left, const pair<double, size_t> &right) {
                return left.first < right.first;
            });

        const size_t family_size = ordered.size();
        if (family_size == 0)
            continue;

        double harmonic = 0.0;
        for (size_t rank = 1; rank <= family_size; rank++)
            harmonic += 1.0 / (double)rank;

        double running_minimum = 1.0;
        for (size_t rank = family_size; rank > 0; rank--) {
            const pair<double, size_t> &entry = ordered[rank - 1];
            double adjusted = entry.first * (double)family_size * harmonic / (double)rank;
            running_minimum = min(running_minimum, min(1.0, adjusted));
            double original_p_value = results[entry.second].p_value;
            if (isfinite(original_p_value) && original_p_value >= 0.0 && original_p_value <= 1.0)
                results[entry.second].fdr_by = running_minimum;
        }

        for (size_t rank = 0; rank < family_size; rank++) {
            SatuTeBranchResult &result = results[ordered[rank].second];
            if (!isfinite(result.p_value) || result.p_value < 0.0 || result.p_value > 1.0)
                continue;
            result.decision_fdr = (result.fdr_by <= fdr_level)
                ? "informative"
                : "saturated";
        }
    }
}

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
    const vector<double> &mode_weights,
    const IntVector *pattern_cat,
    int required_pattern_cat,
    double rate_multiplier,
    int num_threads,
    bool enhanced_weights,
    const SatuTeBranchResult &base,
    const SatuTeDirectionalCache *directional_cache,
    vector<double> *out_pattern_coherence,
    const SatuTeDirectedContext *directed,
    SatuTeMixtureAccumulator *mixture)
{
    return computeCategoryStatisticImpl(
        left_node,
        left_dad,
        right_node,
        right_dad,
        aln,
        model,
        nstates,
        eval,
        formula,
        mode_weights,
        pattern_cat,
        required_pattern_cat,
        rate_multiplier,
        num_threads,
        enhanced_weights,
        base,
        directional_cache,
        out_pattern_coherence,
        directed,
        mixture);
}

} // namespace satute
