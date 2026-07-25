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
    double *partial)
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
            child_partial);

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

    // Rescale to keep deep subtrees away from underflow. Downstream only the
    // normalised posterior partial[s]*freq[s] / sum_s(partial[s]*freq[s]) is
    // used, so any positive scale factor cancels exactly and no log-scale
    // bookkeeping is required. Without this, subtrees of a few hundred taxa
    // drive the product below DBL_MIN and their sites are silently discarded.
    double partial_sum = 0.0;
    for (int state = 0; state < nstates; state++)
        partial_sum += partial[state];
    if (partial_sum > 0.0 && isfinite(partial_sum)) {
        for (int state = 0; state < nstates; state++)
            partial[state] /= partial_sum;
    }
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
    const SatuTeBranchResult &base)
{
    SatuTeBranchResult result = base;
    size_t nmodes = formula.modes.size();
    vector<double> left_second(nmodes * nmodes, 0.0);
    vector<double> right_second(nmodes * nmodes, 0.0);
    vector<double> state_freq(nstates, 0.0);
    // One matrix slot per node id; 2*nseq is a safe upper bound on node count.
    vector<double> transition_cache(
        (2 * aln->getNSeq() + 2) * (size_t)(nstates * nstates), 0.0);
    double *eigenvectors = model->getEigenvectors();

    model->getStateFrequency(&state_freq[0]);

    // Hoisted out of the pattern loop: after this the cache is read-only, which
    // is what makes the loop below safe to share between threads.
    prepareTransitionMatrices(
        left_node, left_dad, model, nstates, rate_multiplier, transition_cache);
    prepareTransitionMatrices(
        right_node, right_dad, model, nstates, rate_multiplier, transition_cache);

    const int nptn = (int)aln->size();

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
    vector<double> thread_coherence(nthreads, 0.0);
    vector<double> thread_valid(nthreads, 0.0);
    vector<int> thread_rate_category_sites(nthreads, 0);
    vector<int> thread_skipped(nthreads, 0);
    vector<double> thread_left_second(nthreads * second_size, 0.0);
    vector<double> thread_right_second(nthreads * second_size, 0.0);

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
    double local_coherence_sum = 0.0;
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
        if (pattern_cat != nullptr && required_pattern_cat >= 0) {
            if ((size_t)ptn >= pattern_cat->size() || (*pattern_cat)[ptn] != required_pattern_cat)
                continue;
        }

        int freq = aln->at(ptn).frequency;
        if (freq <= 0)
            continue;
        local_rate_category_sites += freq;

        computeOrdinaryPartial(
            left_node,
            left_dad,
            ptn,
            aln,
            model,
            nstates,
            transition_cache,
            &scratch[0],
            &left_likelihood[0]);
        computeOrdinaryPartial(
            right_node,
            right_dad,
            ptn,
            aln,
            model,
            nstates,
            transition_cache,
            &scratch[0],
            &right_likelihood[0]);

        double left_sum = 0.0;
        double right_sum = 0.0;
        for (int state = 0; state < nstates; state++) {
            left_posterior[state] = left_likelihood[state] * state_freq[state];
            right_posterior[state] = right_likelihood[state] * state_freq[state];
            left_sum += left_posterior[state];
            right_sum += right_posterior[state];
        }

        if (!isfinite(left_sum) || !isfinite(right_sum) ||
            fabs(left_sum) <= DBL_MIN || fabs(right_sum) <= DBL_MIN) {
            local_skipped_sites += freq;
            continue;
        }
        for (int state = 0; state < nstates; state++) {
            left_posterior[state] /= left_sum;
            right_posterior[state] /= right_sum;
        }

        double site_coherence = 0.0;
        bool valid = true;
        for (size_t i = 0; i < nmodes; i++) {
            left_factor[i] = 0.0;
            right_factor[i] = 0.0;
            for (int state = 0; state < nstates; state++) {
                double eigenvector_value = eigenvectors[state * nstates + formula.modes[i]];
                left_factor[i] += eigenvector_value * left_posterior[state];
                right_factor[i] += eigenvector_value * right_posterior[state];
            }
            if (!isfinite(left_factor[i]) || !isfinite(right_factor[i])) {
                valid = false;
                break;
            }
            site_coherence += mode_weights[i] * left_factor[i] * right_factor[i];
        }
        if (!valid) {
            local_skipped_sites += freq;
            continue;
        }

        local_coherence_sum += freq * site_coherence;
        local_valid_sites += freq;
        for (size_t i = 0; i < nmodes; i++) {
            for (size_t j = 0; j < nmodes; j++) {
                size_t idx = i * nmodes + j;
                local_left_second[idx] += freq * left_factor[i] * left_factor[j];
                local_right_second[idx] += freq * right_factor[i] * right_factor[j];
            }
        }
    }

    thread_coherence[thread_id] = local_coherence_sum;
    thread_valid[thread_id] = local_valid_sites;
    thread_rate_category_sites[thread_id] = local_rate_category_sites;
    thread_skipped[thread_id] = local_skipped_sites;
    for (size_t idx = 0; idx < second_size; idx++) {
        thread_left_second[thread_id * second_size + idx] = local_left_second[idx];
        thread_right_second[thread_id * second_size + idx] = local_right_second[idx];
    }
    } // omp parallel

    // Fixed combination order, so the totals do not depend on thread timing.
    double coherence_sum = 0.0;
    double valid_sites = 0.0;
    int rate_category_sites = 0;
    int skipped_sites = 0;
    for (int thread = 0; thread < nthreads; thread++) {
        coherence_sum += thread_coherence[thread];
        valid_sites += thread_valid[thread];
        rate_category_sites += thread_rate_category_sites[thread];
        skipped_sites += thread_skipped[thread];
        for (size_t idx = 0; idx < second_size; idx++) {
            left_second[idx] += thread_left_second[thread * second_size + idx];
            right_second[idx] += thread_right_second[thread * second_size + idx];
        }
    }

    result.valid_sites = (int)valid_sites;
    result.rate_category_sites = rate_category_sites;
    result.skipped_sites = skipped_sites;
    if (valid_sites <= 0.0) {
        finalizeSatuTeResult(result);
        return result;
    }

    result.coherence = coherence_sum / valid_sites;
    result.variance = 0.0;
    for (size_t i = 0; i < nmodes; i++) {
        for (size_t j = 0; j < nmodes; j++) {
            size_t idx = i * nmodes + j;
            result.variance += mode_weights[i] * mode_weights[j] *
                (left_second[idx] / valid_sites) *
                (right_second[idx] / valid_sites);
        }
    }
    finalizeSatuTeResult(result);
    return result;
}

SatuTeBranchResult poolCategoryResults(
    const vector<SatuTeBranchResult> &category_results,
    const SatuTeBranchResult &base)
{
    SatuTeBranchResult pooled = base;
    pooled.rate_category = "pooled";
    pooled.rate_multiplier = "NA";
    pooled.effective_length = pooled.length;

    double valid_total = 0.0;
    int rate_category_sites = 0;
    int skipped_total = 0;
    for (size_t i = 0; i < category_results.size(); i++) {
        valid_total += category_results[i].valid_sites;
        rate_category_sites += category_results[i].rate_category_sites;
        skipped_total += category_results[i].skipped_sites;
    }

    pooled.valid_sites = (int)valid_total;
    pooled.rate_category_sites = rate_category_sites;
    pooled.skipped_sites = skipped_total;

    if (valid_total <= 0.0) {
        finalizeSatuTeResult(pooled);
        return pooled;
    }

    pooled.coherence = 0.0;
    pooled.variance = 0.0;
    for (size_t i = 0; i < category_results.size(); i++) {
        const SatuTeBranchResult &category = category_results[i];
        if (category.valid_sites <= 0 || !isfinite(category.coherence) || !isfinite(category.variance))
            continue;
        double category_weight = (double)category.valid_sites / valid_total;
        pooled.coherence += category_weight * category.coherence;
        pooled.variance += category_weight * category.variance;
    }
    pooled.weights = "NA";

    finalizeSatuTeResult(pooled);
    return pooled;
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
    const SatuTeBranchResult &base)
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
        base);
}

} // namespace satute
