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

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <limits>
#include <map>
#include <utility>
#include <vector>

using namespace std;

namespace satute {

static void ensureTransitionMatrix(
    Neighbor *branch,
    ModelSubst *model,
    int nstates,
    double rate_multiplier,
    map<Neighbor*, vector<double> > &transition_cache)
{
    if (transition_cache.find(branch) != transition_cache.end())
        return;

    double *eigenvectors = model->getEigenvectors();
    double *inv_eigenvectors = model->getInverseEigenvectors();
    double *eval = model->getEigenvalues();
    vector<double> transition(nstates * nstates, 0.0);
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

    transition_cache[branch] = transition;
}

static void computeOrdinaryPartial(
    Node *node,
    Node *dad,
    size_t ptn,
    Alignment *aln,
    ModelSubst *model,
    int nstates,
    double rate_multiplier,
    map<Neighbor*, vector<double> > &transition_cache,
    vector<double> &partial)
{
    fill(partial.begin(), partial.end(), 0.0);

    if (node->isLeaf()) {
        int state = (node->id >= 0 && node->id < aln->getNSeq()) ? aln->at(ptn)[node->id] : nstates;
        if (state >= 0 && state < nstates) {
            partial[state] = 1.0;
        } else {
            fill(partial.begin(), partial.end(), 1.0);
        }
        return;
    }

    fill(partial.begin(), partial.end(), 1.0);
    FOR_NEIGHBOR_IT(node, dad, it) {
        Neighbor *child_branch = *it;
        vector<double> child_partial(nstates, 0.0);
        computeOrdinaryPartial(
            child_branch->node,
            node,
            ptn,
            aln,
            model,
            nstates,
            rate_multiplier,
            transition_cache,
            child_partial);

        ensureTransitionMatrix(child_branch, model, nstates, rate_multiplier, transition_cache);
        const vector<double> &transition = transition_cache[child_branch];
        for (int parent_state = 0; parent_state < nstates; parent_state++) {
            double contribution = 0.0;
            for (int child_state = 0; child_state < nstates; child_state++) {
                contribution += transition[parent_state * nstates + child_state] *
                    child_partial[child_state];
            }
            partial[parent_state] *= contribution;
        }
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
    const SatuTeBranchResult &base)
{
    SatuTeBranchResult result = base;
    size_t nmodes = formula.modes.size();
    vector<double> left_factor(nmodes, 0.0);
    vector<double> right_factor(nmodes, 0.0);
    vector<double> left_second(nmodes * nmodes, 0.0);
    vector<double> right_second(nmodes * nmodes, 0.0);
    vector<double> state_freq(nstates, 0.0);
    vector<double> left_likelihood(nstates, 0.0);
    vector<double> right_likelihood(nstates, 0.0);
    vector<double> left_posterior(nstates, 0.0);
    vector<double> right_posterior(nstates, 0.0);
    map<Neighbor*, vector<double> > transition_cache;
    double *eigenvectors = model->getEigenvectors();

    model->getStateFrequency(&state_freq[0]);

    double coherence_sum = 0.0;
    double valid_sites = 0.0;
    int rate_category_sites = 0;
    int skipped_sites = 0;

    for (size_t ptn = 0; ptn < aln->size(); ptn++) {
        if (pattern_cat != nullptr && required_pattern_cat >= 0) {
            if (ptn >= pattern_cat->size() || (*pattern_cat)[ptn] != required_pattern_cat)
                continue;
        }

        int freq = aln->at(ptn).frequency;
        if (freq <= 0)
            continue;
        rate_category_sites += freq;

        computeOrdinaryPartial(
            left_node,
            left_dad,
            ptn,
            aln,
            model,
            nstates,
            rate_multiplier,
            transition_cache,
            left_likelihood);
        computeOrdinaryPartial(
            right_node,
            right_dad,
            ptn,
            aln,
            model,
            nstates,
            rate_multiplier,
            transition_cache,
            right_likelihood);

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
            skipped_sites += freq;
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
            skipped_sites += freq;
            continue;
        }

        coherence_sum += freq * site_coherence;
        valid_sites += freq;
        for (size_t i = 0; i < nmodes; i++) {
            for (size_t j = 0; j < nmodes; j++) {
                size_t idx = i * nmodes + j;
                left_second[idx] += freq * left_factor[i] * left_factor[j];
                right_second[idx] += freq * right_factor[i] * right_factor[j];
            }
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
        base);
}

} // namespace satute
