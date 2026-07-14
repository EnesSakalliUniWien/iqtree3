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
#include "satute/satute_report_writer.h"
#include "satute/satute_statistics.h"
#include "satute/satute_support.h"

#include <cstring>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <vector>

using namespace std;
using namespace satute;

void PhyloTree::computeSatuTe(const char *prefix, double alpha, const char *edges_file) {
    if (prefix == nullptr || strlen(prefix) == 0)
        outError("SatuTe requires a valid output prefix");
    if (alpha <= 0.0 || alpha >= 1.0)
        outError("SatuTe alpha must be between 0 and 1");
    if (isSuperTree())
        outError("SatuTe currently supports single-alignment analyses only");
    if (rooted)
        outError("SatuTe currently supports unrooted trees only");
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
    vector<int> all_nonzero_modes = findNonZeroModes(eval, nstates, zero_index);

    DoubleVector pattern_rates;
    IntVector pattern_cat;
    vector<SatuTeRateCategory> rate_categories = buildRateCategories(site_rate, pattern_rates, pattern_cat);
    const IntVector *pattern_cat_ptr = pattern_cat.empty() ? nullptr : &pattern_cat;

    set<int> target_edges = readEdgeFile(edges_file);
    bool restrict_edges = (edges_file != nullptr && strlen(edges_file) > 0);

    BranchVector branches;
    getBranches(branches);
    if (branches.empty())
        outError("SatuTe requires at least one branch");

    cout << "Computing SatuTe branch saturation statistics..." << endl;
    double start_time = getRealTime();

    vector<SatuTeBranchResult> results;
    map<int, Neighbor*> annotated_branches;
    set<int> matched_edges;

    for (BranchVector::iterator brit = branches.begin(); brit != branches.end(); brit++) {
        int branch_id = brit->second->id;
        if (restrict_edges && target_edges.find(branch_id) == target_edges.end())
            continue;

        if (restrict_edges)
            matched_edges.insert(branch_id);

        Neighbor *annotated_branch = brit->second->findNeighbor(brit->first);
        annotated_branches[branch_id] = annotated_branch;

        string label = brit->second->name;
        if (!label.empty())
            annotated_branch->putAttr("label", label);

        int left_taxa = countTipsAwayFrom(brit->first, brit->second);
        int right_taxa = countTipsAwayFrom(brit->second, brit->first);
        string split = splitLabel(brit->first, brit->second);
        double pooled_information_fraction = computePooledInformationFraction(
            eval,
            all_nonzero_modes,
            rate_categories,
            annotated_branch->length);

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
                    formatDouble(category.rate),
                    formatIntVector(formula.modes),
                    eigenvalueString(eval, formula.modes),
                    formatVector(weights),
                    label,
                    split);
                setSaturationScale(
                    base,
                    computeInformationFraction(eval, all_nonzero_modes, effective_length));

                category_results.push_back(
                    computeCategoryStatistic(
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
                formatIntVector(formula.modes),
                eigenvalueString(eval, formula.modes),
                "NA",
                label,
                split);
            setSaturationScale(pooled_base, pooled_information_fraction);

            SatuTeBranchResult pooled = (rate_categories.size() == 1 && rate_categories[0].pattern_cat < 0)
                ? category_results[0]
                : poolCategoryResults(category_results, pooled_base);
            pooled.rate_category = "pooled";
            pooled.rate_multiplier = "NA";
            pooled.effective_length = pooled.length;
            setSaturationScale(pooled, pooled_information_fraction);
            results.push_back(pooled);

            if (!(rate_categories.size() == 1 && rate_categories[0].pattern_cat < 0)) {
                for (size_t category_index = 0; category_index < category_results.size(); category_index++) {
                    if (category_results[category_index].rate_category_sites > 0)
                        results.push_back(category_results[category_index]);
                }
            }
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

    applySeparateFormulaFdr(results, alpha);
    for (size_t result_index = 0; result_index < results.size(); result_index++) {
        const SatuTeBranchResult &result = results[result_index];
        if (result.rate_category != "pooled")
            continue;

        map<int, Neighbor*>::iterator branch_it = annotated_branches.find(result.id);
        if (branch_it == annotated_branches.end())
            continue;
        Neighbor *annotated_branch = branch_it->second;

        if (result.formula == "dominant") {
            annotated_branch->putAttr("satDominantFDR", formatDouble(result.fdr_by));
            annotated_branch->putAttr("satDominantFDRDecision", result.decision_fdr);
            continue;
        }
        if (result.formula != "eigenvalue_weighted")
            continue;

        annotated_branch->putAttr("satFormula", result.formula);
        annotated_branch->putAttr("satC", formatDouble(result.coherence));
        annotated_branch->putAttr("satSE", formatDouble(result.se));
        annotated_branch->putAttr("satZ", formatDouble(result.z_score));
        annotated_branch->putAttr("satP", formatDouble(result.p_value));
        annotated_branch->putAttr("satFDR", formatDouble(result.fdr_by));
        annotated_branch->putAttr("satFDRDecision", result.decision_fdr);
        annotated_branch->putAttr("satInfo", formatDouble(result.information_fraction));
        annotated_branch->putAttr("satIndex", formatDouble(result.saturation_index));
        annotated_branch->putAttr("sat", result.decision);
        annotated_branch->putAttr("satBonf", result.decision_bonferroni);
    }

    ReportWriter().write(*this, prefix, results);
}
