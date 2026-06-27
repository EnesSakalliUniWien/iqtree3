/***************************************************************************
 *   Copyright (C) 2026 by                                                *
 *   IQ-TREE developers                                                    *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 ***************************************************************************/

#include "phylotree.h"
#include "gsl/mygsl.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <vector>

using namespace std;

struct SatuTeFormulaSpec {
    string name;
    vector<int> modes;
    bool eigenvalue_weighted;
};

struct SatuTeRateCategory {
    string label;
    int pattern_cat;
    double rate;
};

struct SatuTeBranchResult {
    int id;
    int left_taxa;
    int right_taxa;
    int valid_sites;
    int skipped_sites;
    double length;
    double effective_length;
    double coherence;
    double variance;
    double se;
    double z_score;
    double p_value;
    double alpha;
    double alpha_adjusted;
    string formula;
    string rate_category;
    string rate_multiplier;
    int rate_category_sites;
    string modes;
    string eigenvalues;
    string weights;
    string decision;
    string decision_bonferroni;
    string label;
    string split;
};

static string satuteDouble(double value) {
    if (!isfinite(value))
        return "NA";
    stringstream ss;
    ss << setprecision(10) << value;
    return ss.str();
}

static string satuteVector(const vector<double> &values) {
    stringstream ss;
    ss << setprecision(10);
    for (size_t i = 0; i < values.size(); i++) {
        if (i > 0)
            ss << ',';
        ss << values[i];
    }
    return ss.str();
}

static string satuteIntVector(const vector<int> &values) {
    stringstream ss;
    for (size_t i = 0; i < values.size(); i++) {
        if (i > 0)
            ss << ',';
        ss << values[i];
    }
    return ss.str();
}

static int countTipsAwayFrom(Node *node, Node *dad) {
    if (node->isLeaf())
        return 1;

    int count = 0;
    FOR_NEIGHBOR_IT(node, dad, it) {
        count += countTipsAwayFrom((*it)->node, node);
    }
    return count;
}

static void collectTipsAwayFrom(Node *node, Node *dad, vector<string> &tips) {
    if (node->isLeaf()) {
        tips.push_back(node->name);
        return;
    }

    FOR_NEIGHBOR_IT(node, dad, it) {
        collectTipsAwayFrom((*it)->node, node, tips);
    }
}

static string satuteSplitLabel(Node *left, Node *right) {
    vector<string> left_tips;
    vector<string> right_tips;
    collectTipsAwayFrom(left, right, left_tips);
    collectTipsAwayFrom(right, left, right_tips);

    vector<string> &chosen = (left_tips.size() <= right_tips.size()) ? left_tips : right_tips;
    sort(chosen.begin(), chosen.end());

    stringstream ss;
    for (size_t i = 0; i < chosen.size(); i++) {
        if (i > 0)
            ss << ',';
        ss << chosen[i];
    }
    return ss.str();
}

static set<int> readSatuTeEdgeFile(const char *filename) {
    set<int> edges;
    if (filename == nullptr || strlen(filename) == 0)
        return edges;

    ifstream in(filename);
    if (!in.is_open())
        outError("Cannot read SatuTe branch ID file ", filename);

    string line;
    while (getline(in, line)) {
        size_t comment = line.find('#');
        if (comment != string::npos)
            line.erase(comment);

        stringstream ss(line);
        string token;
        while (ss >> token) {
            char *end = nullptr;
            long id = strtol(token.c_str(), &end, 10);
            if (end == token.c_str() || *end != '\0')
                outError("Invalid SatuTe branch ID in file ", filename);
            if (id < numeric_limits<int>::min() || id > numeric_limits<int>::max())
                outError("SatuTe branch ID out of integer range in file ", filename);
            edges.insert(id);
        }
    }

    if (edges.empty())
        outError("SatuTe branch ID file is empty");

    return edges;
}

static int findZeroEigenIndex(const double *eval, int nstates) {
    int zero_index = 0;
    double best_abs = fabs(eval[0]);
    for (int i = 1; i < nstates; i++) {
        double value = fabs(eval[i]);
        if (value < best_abs) {
            best_abs = value;
            zero_index = i;
        }
    }
    return zero_index;
}

static vector<int> findNonZeroModes(const double *eval, int nstates, int zero_index) {
    vector<int> modes;
    const double zero_tol = 1e-10;
    for (int i = 0; i < nstates; i++) {
        if (i == zero_index || fabs(eval[i]) <= zero_tol)
            continue;
        modes.push_back(i);
    }
    return modes;
}

static vector<int> findDominantModes(const double *eval, int nstates, int zero_index) {
    vector<int> modes;
    double lambda = -numeric_limits<double>::infinity();
    const double zero_tol = 1e-10;

    for (int i = 0; i < nstates; i++) {
        if (i == zero_index || fabs(eval[i]) <= zero_tol)
            continue;
        if (eval[i] > lambda)
            lambda = eval[i];
    }

    if (!isfinite(lambda))
        return modes;

    double tol = max(1e-8, fabs(lambda) * 1e-6);
    for (int i = 0; i < nstates; i++) {
        if (i == zero_index)
            continue;
        if (fabs(eval[i] - lambda) <= tol)
            modes.push_back(i);
    }
    return modes;
}

static vector<SatuTeFormulaSpec> buildFormulaSpecs(const double *eval, int nstates, int zero_index) {
    vector<int> dominant_modes = findDominantModes(eval, nstates, zero_index);
    vector<int> all_modes = findNonZeroModes(eval, nstates, zero_index);
    if (dominant_modes.empty() || all_modes.empty())
        outError("SatuTe cannot identify non-zero substitution-model eigenvalues");

    vector<SatuTeFormulaSpec> formulas;
    SatuTeFormulaSpec dominant;
    dominant.name = "dominant";
    dominant.modes = dominant_modes;
    dominant.eigenvalue_weighted = false;
    formulas.push_back(dominant);

    SatuTeFormulaSpec weighted;
    weighted.name = "eigenvalue_weighted";
    weighted.modes = all_modes;
    weighted.eigenvalue_weighted = true;
    formulas.push_back(weighted);
    return formulas;
}

static vector<double> buildFormulaWeights(
    const double *eval,
    const SatuTeFormulaSpec &formula,
    double effective_branch_length)
{
    vector<double> weights(formula.modes.size(), 1.0);
    if (!formula.eigenvalue_weighted)
        return weights;

    double dominant = -numeric_limits<double>::infinity();
    for (size_t i = 0; i < formula.modes.size(); i++)
        dominant = max(dominant, eval[formula.modes[i]]);

    for (size_t i = 0; i < formula.modes.size(); i++)
        weights[i] = exp((eval[formula.modes[i]] - dominant) * effective_branch_length);
    return weights;
}

static string eigenvalueString(const double *eval, const vector<int> &modes) {
    vector<double> values;
    for (size_t i = 0; i < modes.size(); i++)
        values.push_back(eval[modes[i]]);
    return satuteVector(values);
}

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

static SatuTeBranchResult makeBaseResult(
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
    result.coherence = numeric_limits<double>::quiet_NaN();
    result.variance = numeric_limits<double>::quiet_NaN();
    result.se = numeric_limits<double>::quiet_NaN();
    result.z_score = numeric_limits<double>::quiet_NaN();
    result.p_value = numeric_limits<double>::quiet_NaN();
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
    result.label = label;
    result.split = split;
    return result;
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

static SatuTeBranchResult computeSatuTeCategoryStatistic(
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

static SatuTeBranchResult poolSatuTeCategoryResults(
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

static vector<SatuTeRateCategory> buildRateCategories(
    RateHeterogeneity *site_rate,
    DoubleVector &pattern_rates,
    IntVector &pattern_cat)
{
    vector<SatuTeRateCategory> categories;
    bool single_rate =
        site_rate->getNRate() == 1 &&
        fabs(site_rate->getRate(0) - 1.0) <= 1e-8 &&
        fabs(site_rate->getProp(0) - 1.0) <= 1e-8 &&
        fabs(site_rate->getPInvar()) <= 1e-12;

    if (single_rate) {
        SatuTeRateCategory category;
        category.label = "pooled";
        category.pattern_cat = -1;
        category.rate = 1.0;
        categories.push_back(category);
        return categories;
    }

    // Keep rate-category assignment identical to IQ-TREE's empirical-Bayes
    // site-rate path: best posterior category per pattern, with IQ-TREE's
    // fitted category multipliers from getRate().
    int observed_categories = site_rate->computePatternRates(pattern_rates, pattern_cat);
    if (pattern_cat.empty()) {
        if (site_rate->getPInvar() > 0.0) {
            outError("SatuTe rate-category analysis requires IQ-TREE pattern-rate assignments. "
                     "IQ-TREE does not emit per-site category rows for a pure +I model; use +I with +G/+R or a discrete rate model.");
        }
        outError("SatuTe could not obtain IQ-TREE pattern-rate assignments for this rate model");
    }

    if (site_rate->getPInvar() > 0.0) {
        SatuTeRateCategory invariant;
        invariant.label = "0";
        invariant.pattern_cat = 0;
        invariant.rate = 0.0;
        categories.push_back(invariant);
        for (int c = 0; c < site_rate->getNRate(); c++) {
            SatuTeRateCategory category;
            category.label = convertIntToString(c + 1);
            category.pattern_cat = c + 1;
            category.rate = site_rate->getRate(c);
            categories.push_back(category);
        }
    } else {
        for (int c = 0; c < site_rate->getNRate(); c++) {
            SatuTeRateCategory category;
            category.label = convertIntToString(c + 1);
            category.pattern_cat = c;
            category.rate = site_rate->getRate(c);
            categories.push_back(category);
        }
    }

    if (observed_categories != (int)categories.size()) {
        outWarning("SatuTe rate-category count differs from IQ-TREE pattern-rate assignment count");
    }
    return categories;
}

void PhyloTree::computeSatuTe(const char *prefix, double alpha, const char *edges_file) {
    if (prefix == nullptr || strlen(prefix) == 0)
        outError("SatuTe requires a valid output prefix");
    if (alpha <= 0.0 || alpha >= 1.0)
        outError("SatuTe alpha must be between 0 and 1");
    if (isSuperTree())
        outError("SatuTe currently supports single-alignment analyses only");
    if (!model || !site_rate || !model_factory)
        outError("SatuTe requires an initialized likelihood model");
    if (!model->useRevKernel())
        outError("SatuTe currently requires the reversible likelihood kernel");
    if (model->getNMixtures() != 1)
        outError("SatuTe currently does not support mixture models");
    if (model_factory->getASC() != ASC_NONE || !model_factory->unobserved_ptns.empty())
        outError("SatuTe currently does not support ascertainment bias correction");

    int nstates = aln->num_states;
    if (nstates < 2)
        outError("SatuTe requires at least two model states");

    double *eval = model->getEigenvalues();
    if (eval == nullptr)
        outError("SatuTe requires eigenvalues from the substitution model");

    int zero_index = findZeroEigenIndex(eval, nstates);
    vector<SatuTeFormulaSpec> formulas = buildFormulaSpecs(eval, nstates, zero_index);

    DoubleVector pattern_rates;
    IntVector pattern_cat;
    vector<SatuTeRateCategory> rate_categories = buildRateCategories(site_rate, pattern_rates, pattern_cat);
    const IntVector *pattern_cat_ptr = pattern_cat.empty() ? nullptr : &pattern_cat;

    set<int> target_edges = readSatuTeEdgeFile(edges_file);
    bool restrict_edges = (edges_file != nullptr && strlen(edges_file) > 0);

    BranchVector branches;
    getBranches(branches);
    if (branches.empty())
        outError("SatuTe requires at least one branch");

    cout << "Computing SatuTe branch saturation statistics..." << endl;
    double start_time = getRealTime();

    vector<SatuTeBranchResult> results;
    set<int> matched_edges;

    for (BranchVector::iterator brit = branches.begin(); brit != branches.end(); brit++) {
        int branch_id = brit->second->id;
        if (restrict_edges && target_edges.find(branch_id) == target_edges.end())
            continue;

        if (restrict_edges)
            matched_edges.insert(branch_id);

        Neighbor *annotated_branch = brit->second->findNeighbor(brit->first);

        string label = brit->second->name;
        if (!label.empty())
            annotated_branch->putAttr("label", label);

        int left_taxa = countTipsAwayFrom(brit->first, brit->second);
        int right_taxa = countTipsAwayFrom(brit->second, brit->first);
        string split = satuteSplitLabel(brit->first, brit->second);
        SatuTeBranchResult dominant_pooled_for_tree;
        bool have_dominant_pooled = false;

        for (size_t formula_index = 0; formula_index < formulas.size(); formula_index++) {
            const SatuTeFormulaSpec &formula = formulas[formula_index];
            vector<SatuTeBranchResult> category_results;

            for (size_t category_index = 0; category_index < rate_categories.size(); category_index++) {
                const SatuTeRateCategory &category = rate_categories[category_index];
                double effective_length = annotated_branch->length * category.rate;
                vector<double> weights = buildFormulaWeights(eval, formula, effective_length);
                SatuTeBranchResult base = makeBaseResult(
                    branch_id,
                    left_taxa,
                    right_taxa,
                    annotated_branch->length,
                    effective_length,
                    alpha,
                    formula.name,
                    category.label,
                    satuteDouble(category.rate),
                    satuteIntVector(formula.modes),
                    eigenvalueString(eval, formula.modes),
                    satuteVector(weights),
                    label,
                    split);

                category_results.push_back(
                    computeSatuTeCategoryStatistic(
                        brit->first,
                        brit->second,
                        brit->second,
                        brit->first,
                        aln,
                        model,
                        nstates,
                        eval,
                        formula,
                        weights,
                        pattern_cat_ptr,
                        category.pattern_cat,
                        category.rate,
                        base));
            }

            SatuTeBranchResult pooled_base = makeBaseResult(
                branch_id,
                left_taxa,
                right_taxa,
                annotated_branch->length,
                annotated_branch->length,
                alpha,
                formula.name,
                "pooled",
                "NA",
                satuteIntVector(formula.modes),
                eigenvalueString(eval, formula.modes),
                "NA",
                label,
                split);

            SatuTeBranchResult pooled = (rate_categories.size() == 1 && rate_categories[0].pattern_cat < 0)
                ? category_results[0]
                : poolSatuTeCategoryResults(category_results, pooled_base);
            pooled.rate_category = "pooled";
            pooled.rate_multiplier = "NA";
            pooled.effective_length = pooled.length;
            results.push_back(pooled);

            if (formula.name == "dominant") {
                dominant_pooled_for_tree = pooled;
                have_dominant_pooled = true;
            }

            if (!(rate_categories.size() == 1 && rate_categories[0].pattern_cat < 0)) {
                for (size_t category_index = 0; category_index < category_results.size(); category_index++)
                    results.push_back(category_results[category_index]);
            }
        }

        if (have_dominant_pooled) {
            annotated_branch->putAttr("satC", satuteDouble(dominant_pooled_for_tree.coherence));
            annotated_branch->putAttr("satSE", satuteDouble(dominant_pooled_for_tree.se));
            annotated_branch->putAttr("satZ", satuteDouble(dominant_pooled_for_tree.z_score));
            annotated_branch->putAttr("satP", satuteDouble(dominant_pooled_for_tree.p_value));
            annotated_branch->putAttr("sat", dominant_pooled_for_tree.decision);
            annotated_branch->putAttr("satBonf", dominant_pooled_for_tree.decision_bonferroni);
        }
    }

    cout << getRealTime() - start_time << " sec" << endl;
    if (restrict_edges) {
        vector<int> missing_edges;
        for (set<int>::iterator it = target_edges.begin(); it != target_edges.end(); it++) {
            if (matched_edges.find(*it) == matched_edges.end())
                missing_edges.push_back(*it);
        }
        if (!missing_edges.empty()) {
            stringstream ss;
            ss << "Requested SatuTe branch IDs not found:";
            for (size_t i = 0; i < missing_edges.size(); i++)
                ss << ' ' << missing_edges[i];
            outError(ss.str());
        }
    }

    string out_prefix = prefix;
    string filename = out_prefix + ".sat.stat";
    string str = out_prefix + ".sat.tree";
    printTree(str.c_str(), WT_BR_LEN + WT_NEWLINE);
    cout << "Tree with SatuTe annotations written to " << str << endl;

    str = out_prefix + ".sat.tree.nex";
    printNexus(str, WT_BR_LEN, "See " + filename + " for SatuTe branch saturation statistics. "
               "This file is best viewed in FigTree.");
    cout << "Annotated tree (best viewed in FigTree) written to " << str << endl;

    str = out_prefix + ".sat.branch";
    printTree(str.c_str(), WT_BR_LEN + WT_INT_NODE + WT_NEWLINE);
    cout << "Tree with SatuTe branch IDs written to " << str << endl;

    ofstream out(filename.c_str());
    out << "# SatuTe branch saturation statistics" << endl
        << "# This file can be read in MS Excel or in R with command:" << endl
        << "#   tab=read.table('" << filename << "',header=TRUE)" << endl
        << "# Columns are tab-separated with following meaning:" << endl
        << "#   ID: Branch ID" << endl
        << "#   Formula: dominant or eigenvalue_weighted" << endl
        << "#   RateCategory: pooled for all sites, or the IQ-TREE maximum-posterior pattern-rate category" << endl
        << "#   RateMultiplier: IQ-TREE category rate used to rescale branch lengths; NA for pooled rows" << endl
        << "#   RateSites: Number of sites assigned by IQ-TREE to this rate category or pooled row" << endl
        << "#   LeftTaxa, RightTaxa: Number of taxa on each side of the branch" << endl
        << "#   ValidSites: Alignment sites used after skipping invalid likelihood ratios" << endl
        << "#   SkippedSites: Alignment sites skipped because likelihood ratios were invalid" << endl
        << "#   satC: Mean SatuTe branch coherence statistic" << endl
        << "#   satVar: Null variance estimate of satC" << endl
        << "#   satSE: Standard error of satC" << endl
        << "#   satZ: One-sided normal test statistic" << endl
        << "#   satP: One-sided p-value; low values indicate phylogenetic signal across the branch" << endl
        << "#   Alpha: Unadjusted significance level" << endl
        << "#   AlphaTaxonBonf: Alpha adjusted by LeftTaxa*RightTaxa" << endl
        << "#   Decision: informative if satP <= Alpha, otherwise saturated" << endl
        << "#   DecisionTaxonBonf: informative if satP <= AlphaTaxonBonf" << endl
        << "#   Label: Existing branch label" << endl
        << "#   Length: Original branch length" << endl
        << "#   EffectiveLength: Branch length after category-rate rescaling" << endl
        << "#   Split: Comma-separated taxa on the smaller side of the branch" << endl
        << "#   Modes, Eigenvalues, Weights: spectral modes used by the formula; pooled rate-category rows use NA weights because category-specific weights are listed in category rows" << endl
        << "ID\tFormula\tRateCategory\tRateMultiplier\tRateSites\tLeftTaxa\tRightTaxa\tValidSites\tSkippedSites\tsatC\tsatVar\tsatSE\tsatZ\tsatP\tAlpha\tAlphaTaxonBonf\tDecision\tDecisionTaxonBonf\tLabel\tLength\tEffectiveLength\tSplit\tModes\tEigenvalues\tWeights" << endl;

    for (size_t i = 0; i < results.size(); i++) {
        SatuTeBranchResult &result = results[i];
        out << result.id
            << '\t' << result.formula
            << '\t' << result.rate_category
            << '\t' << result.rate_multiplier
            << '\t' << result.rate_category_sites
            << '\t' << result.left_taxa
            << '\t' << result.right_taxa
            << '\t' << result.valid_sites
            << '\t' << result.skipped_sites
            << '\t' << satuteDouble(result.coherence)
            << '\t' << satuteDouble(result.variance)
            << '\t' << satuteDouble(result.se)
            << '\t' << satuteDouble(result.z_score)
            << '\t' << satuteDouble(result.p_value)
            << '\t' << satuteDouble(result.alpha)
            << '\t' << satuteDouble(result.alpha_adjusted)
            << '\t' << result.decision
            << '\t' << result.decision_bonferroni
            << '\t' << result.label
            << '\t' << satuteDouble(result.length)
            << '\t' << satuteDouble(result.effective_length)
            << '\t' << result.split
            << '\t' << result.modes
            << '\t' << result.eigenvalues
            << '\t' << result.weights
            << endl;
    }
    out.close();
    cout << "SatuTe statistics per branch printed to " << filename << endl;
}
