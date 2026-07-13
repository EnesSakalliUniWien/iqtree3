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

static void computeScaledPartial(
    Node *node,
    Node *dad,
    size_t ptn,
    Alignment *aln,
    ModelSubst *model,
    int nstates,
    double rate_multiplier,
    map<Neighbor*, vector<double> > &transition_cache,
    vector<double> &partial,
    double &log_scale)
{
    fill(partial.begin(), partial.end(), 0.0);
    log_scale = 0.0;

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
        double child_log_scale = 0.0;
        computeScaledPartial(
            child_branch->node,
            node,
            ptn,
            aln,
            model,
            nstates,
            rate_multiplier,
            transition_cache,
            child_partial,
            child_log_scale);

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
        log_scale += child_log_scale;

        double scale = 0.0;
        for (int state = 0; state < nstates; state++)
            scale = max(scale, fabs(partial[state]));
        if (scale > 0.0 && isfinite(scale)) {
            for (int state = 0; state < nstates; state++)
                partial[state] /= scale;
            log_scale += log(scale);
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

static SatuTeBranchResult computeMixtureStatisticImpl(
    Node *left_node,
    Node *left_dad,
    Node *right_node,
    Node *right_dad,
    Alignment *aln,
    ModelSubst *model,
    int nstates,
    const double *eval,
    const vector<int> &modes,
    const vector<SatuTeRateCategory> &categories,
    double branch_length,
    double log_weight_shift,
    const SatuTeBranchResult &base)
{
    SatuTeBranchResult result = base;
    vector<double> state_freq(nstates, 0.0);
    vector<double> left_partial(nstates, 0.0);
    vector<double> right_partial(nstates, 0.0);
    vector<double> left_posterior(nstates, 0.0);
    vector<double> right_posterior(nstates, 0.0);
    vector<map<Neighbor*, vector<double> > > transition_caches(categories.size());
    double *eigenvectors = model->getEigenvectors();
    model->getStateFrequency(&state_freq[0]);

    double statistic_sum = 0.0;
    double statistic_square_sum = 0.0;
    double valid_sites = 0.0;
    int assigned_sites = 0;
    int skipped_sites = 0;

    for (size_t ptn = 0; ptn < aln->size(); ptn++) {
        int freq = aln->at(ptn).frequency;
        if (freq <= 0)
            continue;
        assigned_sites += freq;

        vector<double> log_bases;
        vector<double> signals;
        log_bases.reserve(categories.size());
        signals.reserve(categories.size());
        bool site_valid = true;

        for (size_t category_index = 0; category_index < categories.size(); category_index++) {
            const SatuTeRateCategory &category = categories[category_index];
            if (!(category.proportion > 0.0) || !isfinite(category.proportion))
                continue;
            if (category.rate < 0.0 || !isfinite(category.rate)) {
                site_valid = false;
                break;
            }

            double left_log_scale = 0.0;
            double right_log_scale = 0.0;
            computeScaledPartial(
                left_node,
                left_dad,
                ptn,
                aln,
                model,
                nstates,
                category.rate,
                transition_caches[category_index],
                left_partial,
                left_log_scale);
            computeScaledPartial(
                right_node,
                right_dad,
                ptn,
                aln,
                model,
                nstates,
                category.rate,
                transition_caches[category_index],
                right_partial,
                right_log_scale);

            double left_sum = 0.0;
            double right_sum = 0.0;
            for (int state = 0; state < nstates; state++) {
                left_posterior[state] = left_partial[state] * state_freq[state];
                right_posterior[state] = right_partial[state] * state_freq[state];
                left_sum += left_posterior[state];
                right_sum += right_posterior[state];
            }

            if (!isfinite(left_sum) || !isfinite(right_sum) ||
                !isfinite(left_log_scale) || !isfinite(right_log_scale)) {
                site_valid = false;
                break;
            }

            if (!(category.rate > 0.0)) {
                double invariant_joint = 0.0;
                for (int state = 0; state < nstates; state++)
                    invariant_joint += state_freq[state] * left_partial[state] * right_partial[state];
                if (!isfinite(invariant_joint) || invariant_joint < 0.0) {
                    site_valid = false;
                    break;
                }
                if (!(invariant_joint > 0.0))
                    continue;
                log_bases.push_back(
                    log(category.proportion) + log(invariant_joint) +
                    left_log_scale + right_log_scale);
                signals.push_back(0.0);
                continue;
            }

            if (!(left_sum > 0.0) || !(right_sum > 0.0))
                continue;
            for (int state = 0; state < nstates; state++) {
                left_posterior[state] /= left_sum;
                right_posterior[state] /= right_sum;
            }

            double category_signal = 0.0;
            bool category_valid = true;
            double effective_length = branch_length * category.rate;
            for (size_t mode_index = 0; mode_index < modes.size(); mode_index++) {
                double left_factor = 0.0;
                double right_factor = 0.0;
                for (int state = 0; state < nstates; state++) {
                    double eigenvector_value = eigenvectors[state * nstates + modes[mode_index]];
                    left_factor += eigenvector_value * left_posterior[state];
                    right_factor += eigenvector_value * right_posterior[state];
                }
                double globally_scaled_weight = exp(
                    eval[modes[mode_index]] * effective_length - log_weight_shift);
                double contribution = globally_scaled_weight * left_factor * right_factor;
                if (!isfinite(contribution)) {
                    category_valid = false;
                    break;
                }
                category_signal += contribution;
            }
            if (!category_valid) {
                site_valid = false;
                break;
            }

            log_bases.push_back(
                log(category.proportion) +
                log(left_sum) + left_log_scale +
                log(right_sum) + right_log_scale);
            signals.push_back(category_signal);
        }

        if (!site_valid || log_bases.empty()) {
            skipped_sites += freq;
            continue;
        }

        double max_log_base = *max_element(log_bases.begin(), log_bases.end());
        double denominator = 0.0;
        double numerator = 0.0;
        for (size_t category_index = 0; category_index < log_bases.size(); category_index++) {
            double responsibility_numerator = exp(log_bases[category_index] - max_log_base);
            denominator += responsibility_numerator;
            numerator += responsibility_numerator * signals[category_index];
        }
        if (!(denominator > 0.0) || !isfinite(denominator) || !isfinite(numerator)) {
            skipped_sites += freq;
            continue;
        }

        double site_statistic = numerator / denominator;
        if (!isfinite(site_statistic)) {
            skipped_sites += freq;
            continue;
        }
        statistic_sum += freq * site_statistic;
        statistic_square_sum += freq * site_statistic * site_statistic;
        valid_sites += freq;
    }

    result.valid_sites = (int)valid_sites;
    result.rate_category_sites = assigned_sites;
    result.skipped_sites = skipped_sites;
    if (valid_sites <= 1.0) {
        finalizeSatuTeResult(result);
        return result;
    }

    result.coherence = statistic_sum / valid_sites;
    double centered_sum = statistic_square_sum - valid_sites * result.coherence * result.coherence;
    if (centered_sum < 0.0 && centered_sum > -1e-12 * max(1.0, statistic_square_sum))
        centered_sum = 0.0;
    result.variance = centered_sum / (valid_sites - 1.0);
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

namespace {

class CategoryStatisticStrategy : public StatisticStrategy {
public:
    SatuTeBranchResult compute(const StatisticRequest &request) const override {
        return computeCategoryStatisticImpl(
            request.left_node, request.left_dad, request.right_node, request.right_dad,
            request.aln, request.model, request.nstates, request.eigenvalues,
            *request.formula, *request.mode_weights, request.pattern_categories,
            request.required_pattern_category, request.rate_multiplier, *request.base);
    }
};

class MixtureStatisticStrategy : public StatisticStrategy {
public:
    SatuTeBranchResult compute(const StatisticRequest &request) const override {
        return computeMixtureStatisticImpl(
            request.left_node, request.left_dad, request.right_node, request.right_dad,
            request.aln, request.model, request.nstates, request.eigenvalues,
            *request.modes, *request.rate_categories, request.branch_length,
            request.log_weight_shift, *request.base);
    }
};

} // namespace

unique_ptr<StatisticStrategy> StatisticStrategyFactory::create(StatisticKind kind) {
    if (kind == StatisticKind::Category)
        return unique_ptr<StatisticStrategy>(new CategoryStatisticStrategy());
    return unique_ptr<StatisticStrategy>(new MixtureStatisticStrategy());
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
    StatisticRequest request = {
        left_node, left_dad, right_node, right_dad, aln, model, nstates, eval,
        &formula, &mode_weights, pattern_cat, required_pattern_cat, rate_multiplier,
        nullptr, nullptr, 0.0, 0.0, &base
    };
    return StatisticStrategyFactory::create(StatisticKind::Category)->compute(request);
}

SatuTeBranchResult computeMixtureStatistic(
    Node *left_node,
    Node *left_dad,
    Node *right_node,
    Node *right_dad,
    Alignment *aln,
    ModelSubst *model,
    int nstates,
    const double *eval,
    const vector<int> &modes,
    const vector<SatuTeRateCategory> &categories,
    double branch_length,
    double log_weight_shift,
    const SatuTeBranchResult &base)
{
    StatisticRequest request = {
        left_node, left_dad, right_node, right_dad, aln, model, nstates, eval,
        nullptr, nullptr, nullptr, -1, 0.0, &modes, &categories,
        branch_length, log_weight_shift, &base
    };
    return StatisticStrategyFactory::create(StatisticKind::Mixture)->compute(request);
}

} // namespace satute
