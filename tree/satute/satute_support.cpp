/***************************************************************************
 *   Copyright (C) 2026 by                                                *
 *   IQ-TREE developers                                                    *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 ***************************************************************************/

#include "satute_support.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

using namespace std;

namespace satute {

string formatDouble(double value) {
    if (!isfinite(value))
        return "NA";
    stringstream ss;
    ss << setprecision(10) << value;
    return ss.str();
}

string formatVector(const vector<double> &values) {
    stringstream ss;
    ss << setprecision(10);
    for (size_t i = 0; i < values.size(); i++) {
        if (i > 0)
            ss << ',';
        ss << values[i];
    }
    return ss.str();
}

string formatIntVector(const vector<int> &values) {
    stringstream ss;
    for (size_t i = 0; i < values.size(); i++) {
        if (i > 0)
            ss << ',';
        ss << values[i];
    }
    return ss.str();
}

int countTipsAwayFrom(Node *node, Node *dad) {
    if (node->isLeaf())
        return 1;

    int count = 0;
    FOR_NEIGHBOR_IT(node, dad, it) {
        count += countTipsAwayFrom((*it)->node, node);
    }
    return count;
}

void collectTipsAwayFrom(Node *node, Node *dad, vector<string> &tips) {
    if (node->isLeaf()) {
        tips.push_back(node->name);
        return;
    }

    FOR_NEIGHBOR_IT(node, dad, it) {
        collectTipsAwayFrom((*it)->node, node, tips);
    }
}

string splitLabel(Node *left, Node *right) {
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

set<int> readEdgeFile(const char *filename) {
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

int findZeroEigenIndex(const double *eval, int nstates) {
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

vector<int> findNonZeroModes(const double *eval, int nstates, int zero_index) {
    vector<int> modes;
    const double zero_tol = 1e-10;
    for (int i = 0; i < nstates; i++) {
        if (i == zero_index || fabs(eval[i]) <= zero_tol)
            continue;
        modes.push_back(i);
    }
    return modes;
}

vector<int> findDominantModes(const double *eval, int nstates, int zero_index) {
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

vector<SatuTeFormulaSpec> buildFormulaSpecs(const double *eval, int nstates, int zero_index) {
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

vector<double> buildFormulaWeights(
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

string eigenvalueString(const double *eval, const vector<int> &modes) {
    vector<double> values;
    for (size_t i = 0; i < modes.size(); i++)
        values.push_back(eval[modes[i]]);
    return formatVector(values);
}

double computeInformationFraction(
    const double *eval,
    const vector<int> &modes,
    double effective_length)
{
    if (modes.empty() || effective_length < 0.0 || !isfinite(effective_length))
        return numeric_limits<double>::quiet_NaN();

    double information = 0.0;
    for (size_t i = 0; i < modes.size(); i++)
        information += exp(2.0 * eval[modes[i]] * effective_length);
    information /= modes.size();
    if (information < 0.0 && information > -1e-12)
        information = 0.0;
    if (information > 1.0 && information < 1.0 + 1e-12)
        information = 1.0;
    return information;
}

double computeMixtureInformationFraction(
    const double *eval,
    const vector<int> &modes,
    const vector<SatuTeRateCategory> &categories,
    double branch_length)
{
    double weighted_information = 0.0;
    double total_proportion = 0.0;
    for (size_t i = 0; i < categories.size(); i++) {
        const SatuTeRateCategory &category = categories[i];
        if (!(category.proportion > 0.0) || !isfinite(category.proportion))
            continue;
        double category_information = computeInformationFraction(
            eval,
            modes,
            branch_length * category.rate);
        if (!isfinite(category_information))
            return numeric_limits<double>::quiet_NaN();
        weighted_information += category.proportion * category_information;
        total_proportion += category.proportion;
    }
    if (!(total_proportion > 0.0))
        return numeric_limits<double>::quiet_NaN();
    return weighted_information / total_proportion;
}

double computeMixtureLogWeightShift(
    const double *eval,
    const vector<int> &modes,
    const vector<SatuTeRateCategory> &categories,
    double branch_length)
{
    double shift = -numeric_limits<double>::infinity();
    for (size_t category_index = 0; category_index < categories.size(); category_index++) {
        const SatuTeRateCategory &category = categories[category_index];
        if (!(category.proportion > 0.0) || !(category.rate > 0.0) ||
            !isfinite(category.proportion) || !isfinite(category.rate))
            continue;
        double effective_length = branch_length * category.rate;
        for (size_t mode_index = 0; mode_index < modes.size(); mode_index++)
            shift = max(shift, eval[modes[mode_index]] * effective_length);
    }
    return isfinite(shift) ? shift : 0.0;
}

void setSaturationScale(SatuTeBranchResult &result, double information_fraction) {
    result.information_fraction = information_fraction;
    result.saturation_index = isfinite(information_fraction)
        ? 1.0 - information_fraction
        : numeric_limits<double>::quiet_NaN();
    if (result.saturation_index < 0.0 && result.saturation_index > -1e-12)
        result.saturation_index = 0.0;
    if (result.saturation_index > 1.0 && result.saturation_index < 1.0 + 1e-12)
        result.saturation_index = 1.0;
}

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

vector<SatuTeRateCategory> buildRateCategories(
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
        category.proportion = 1.0;
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
        invariant.proportion = site_rate->getPInvar();
        categories.push_back(invariant);
        for (int c = 0; c < site_rate->getNRate(); c++) {
            SatuTeRateCategory category;
            category.label = convertIntToString(c + 1);
            category.pattern_cat = c + 1;
            category.rate = site_rate->getRate(c);
            category.proportion = site_rate->getProp(c);
            categories.push_back(category);
        }
    } else {
        for (int c = 0; c < site_rate->getNRate(); c++) {
            SatuTeRateCategory category;
            category.label = convertIntToString(c + 1);
            category.pattern_cat = c;
            category.rate = site_rate->getRate(c);
            category.proportion = site_rate->getProp(c);
            categories.push_back(category);
        }
    }

    if (observed_categories != (int)categories.size()) {
        outWarning("SatuTe rate-category count differs from IQ-TREE pattern-rate assignment count");
    }
    return categories;
}


} // namespace satute
