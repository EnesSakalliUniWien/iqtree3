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
#include "satute/satute_edges.h"
#include "satute/satute_format.h"
#include "satute/satute_formulas.h"
#include "satute/satute_information.h"
#include "satute/satute_rates.h"
#include "satute/satute_results.h"
#include "satute/satute_topology.h"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <utility>
#include <vector>

using namespace std;
using namespace satute;

void PhyloTree::computeSatuTe(const char *prefix, double alpha, const char *edges_file, bool enhanced_weights, bool write_sites, int sites_max) {
    if (prefix == nullptr || strlen(prefix) == 0)
        outError("SatuTe requires a valid output prefix");
    if (alpha <= 0.0 || alpha >= 1.0)
        outError("SatuTe alpha must be between 0 and 1");
    if (isSuperTree())
        outError("SatuTe currently supports single-alignment analyses only");
    if (!model || !site_rate || !model_factory)
        outError("SatuTe requires an initialized likelihood model");
    // Reversible models keep the published unrooted formulation. Without
    // detailed balance the pulley principle fails, the root position is
    // identifiable and the statistic lives on rooted, directed edges, so the
    // non-reversible path requires the rooted tree IQ-TREE itself maintains
    // for these models.
    const bool directed_nonreversible = !model->isReversible();
    if (directed_nonreversible) {
        if (!rooted)
            outError("SatuTe requires a rooted tree for non-reversible models");
    } else {
        if (rooted)
            outError("SatuTe currently supports unrooted trees only");
        if (!model->useRevKernel())
            outError("SatuTe currently requires the reversible likelihood kernel");
    }
    if (model->getNMixtures() != 1)
        outError("SatuTe currently does not support mixture models");
    if (model_factory->getASC() != ASC_NONE || !model_factory->unobserved_ptns.empty())
        outError("SatuTe currently does not support ascertainment bias correction");

    int nstates = aln->num_states;
    if (nstates < 2)
        outError("SatuTe requires at least two model states");

    double *eval = nullptr;
    vector<SatuTeFormulaSpec> formula_specs;
    vector<int> all_nonzero_modes;
    if (directed_nonreversible) {
        // Mode sets depend on the branch length through the singular basis of
        // M(t), so only the two formula labels are fixed model-wide.
        formula_specs = formulas::buildDirectedFormulaSpecs();
    } else {
        eval = model->getEigenvalues();
        if (eval == nullptr)
            outError("SatuTe requires eigenvalues from the substitution model");
        int zero_index = formulas::findZeroEigenIndex(eval, nstates);
        formula_specs = formulas::buildFormulaSpecs(eval, nstates, zero_index);
        all_nonzero_modes = formulas::findNonZeroModes(eval, nstates, zero_index);
    }

    DoubleVector pattern_rates;
    IntVector pattern_cat;
    vector<SatuTeRateCategory> rate_categories = rates::buildRateCategories(site_rate, pattern_rates, pattern_cat);
    const IntVector *pattern_cat_ptr = pattern_cat.empty() ? nullptr : &pattern_cat;

    set<int> target_edges = edges::readEdgeFile(edges_file);
    bool restrict_edges = (edges_file != nullptr && strlen(edges_file) > 0);

    BranchVector branches;
    getBranches(branches);
    if (branches.empty())
        outError("SatuTe requires at least one branch");

    cout << "Computing SatuTe branch saturation statistics..." << endl;
    double start_time = getRealTime();

    // Every branch needs the two directional partial likelihoods that meet on
    // it. Computed per branch that is O(n^2) in taxa, because each of the 2n-3
    // branches re-walks a whole subtree. Computed once per rate category by a
    // two-pass sweep it is O(n), so build that shared cache when it fits in
    // memory and fall back to the per-branch recursion when it does not.
    const size_t node_bound = 2 * aln->getNSeq() + 2;
    // Parent relation for the directed statistic: which end of every edge is
    // the true parent under the rooted tree.
    vector<int> parent_of;
    if (directed_nonreversible) {
        parent_of.assign(node_bound, -1);
        topology::fillParentMap(root, nullptr, parent_of);
    }
    // Every node is an endpoint of some branch, so this sees the whole tree.
    int max_degree = 0;
    for (BranchVector::iterator brit = branches.begin(); brit != branches.end(); brit++) {
        max_degree = max(max_degree, (int)brit->first->neighbors.size());
        max_degree = max(max_degree, (int)brit->second->neighbors.size());
    }
    Node *sweep_root = nullptr;
    for (BranchVector::iterator brit = branches.begin();
         brit != branches.end() && sweep_root == nullptr; brit++) {
        if (!brit->first->isLeaf())
            sweep_root = brit->first;
        else if (!brit->second->isLeaf())
            sweep_root = brit->second;
    }

    // Rate category is the outer loop below, so only one category's cache is
    // ever resident and the budget covers a single category.
    const size_t directional_cache_budget = (size_t)2 * 1024 * 1024 * 1024;
    const size_t directional_cache_bytes = (max_degree > 0)
        ? satute::directionalCacheBytes(node_bound, max_degree, aln->size(), nstates)
        : 0;
    // Escape hatch: forces the per-branch recursion so the two paths can be
    // compared against each other within one binary. Validation only.
    const bool directional_cache_disabled =
        (getenv("SATUTE_DISABLE_DIRECTIONAL_CACHE") != nullptr);
    // The shared message cache transports every edge with the reversible
    // kernel convention, so the direction-aware non-reversible path always
    // uses the per-branch recursion.
    const bool use_directional_cache =
        !directed_nonreversible &&
        !directional_cache_disabled && (sweep_root != nullptr) &&
        (directional_cache_bytes > 0) &&
        (directional_cache_bytes <= directional_cache_budget);
    if (!directed_nonreversible &&
        !use_directional_cache && sweep_root != nullptr && !directional_cache_disabled) {
        cout << "  Directional partial-likelihood cache would need "
             << (directional_cache_bytes / (1024 * 1024))
             << " MB per rate category; using the slower per-branch recursion"
             << " instead." << endl;
    }

    vector<SatuTeBranchResult> results;
    // Per-site rows, collected only under --satute-sites. Pooled across rate
    // categories, one entry per (branch, formula).
    vector<SatuTeSiteSeries> site_series;
    map<int, Neighbor*> annotated_branches;
    set<int> matched_edges;

    // Everything about a branch that does not depend on the rate category.
    // Resolved once so the compute pass below can run category-outer.
    struct AnalysedBranch {
        Node *left_node;
        Node *right_node;
        Neighbor *annotated_branch;
        int branch_id;
        int left_taxa;
        int right_taxa;
        string label;
        string split;
        double sibling_length;
        double pooled_information_fraction;
    };
    vector<AnalysedBranch> analysed;

    for (BranchVector::iterator brit = branches.begin(); brit != branches.end(); brit++) {
        int branch_id = brit->second->id;
        // The edge to the __root__ leaf splits {} from all taxa: nothing lies
        // above it, so no saturation test is defined there.
        if (directed_nonreversible &&
            (isRootLeaf(brit->first) || isRootLeaf(brit->second)))
            continue;
        if (restrict_edges && target_edges.find(branch_id) == target_edges.end())
            continue;

        if (restrict_edges)
            matched_edges.insert(branch_id);

        Neighbor *annotated_branch = brit->second->findNeighbor(brit->first);
        annotated_branches[branch_id] = annotated_branch;

        string label = brit->second->name;
        if (!label.empty())
            annotated_branch->putAttr("label", label);

        AnalysedBranch entry;
        entry.left_node = brit->first;
        entry.right_node = brit->second;
        // The directed statistic distinguishes the rootward (left) and
        // tipward (right) side; orient without touching branch ids.
        if (directed_nonreversible &&
            parent_of[brit->first->id] == brit->second->id)
            swap(entry.left_node, entry.right_node);
        entry.annotated_branch = annotated_branch;
        entry.branch_id = branch_id;
        entry.left_taxa = topology::countTipsAwayFrom(entry.left_node, entry.right_node);
        entry.right_taxa = topology::countTipsAwayFrom(entry.right_node, entry.left_node);
        entry.label = label;
        entry.split = topology::splitLabel(entry.left_node, entry.right_node);
        // Below a two-child root the rootward side is only the sibling clade,
        // so the tested split spans both root edges.
        entry.sibling_length = 0.0;
        if (directed_nonreversible && parent_of[entry.left_node->id] == root->id &&
            entry.left_node->neighbors.size() == 3) {
            FOR_NEIGHBOR_IT(entry.left_node, entry.right_node, it)
                if (!isRootLeaf((*it)->node))
                    entry.sibling_length = (*it)->length;
        }
        entry.pooled_information_fraction = directed_nonreversible
            ? information::computePooledDirectedInformationFraction(
                  model, nstates, rate_categories, annotated_branch->length,
                  entry.sibling_length)
            : information::computePooledInformationFraction(
                  eval, all_nonzero_modes, rate_categories, annotated_branch->length);
        analysed.push_back(entry);
    }

    // Compute pass. Rate category is outermost so each category's directional
    // cache is built, used by every branch and formula, then released before
    // the next one -- only one is resident at a time. Results are held in a
    // [branch][formula][category] cube and reduced afterwards, which keeps the
    // emitted row order identical to the branch-major original.
    vector<vector<vector<SatuTeBranchResult> > > category_stats(
        analysed.size(),
        vector<vector<SatuTeBranchResult> >(formula_specs.size()));
    // Rate categories partition the patterns, so summing the per-category
    // arrays reconstructs the pooled per-pattern value.
    vector<vector<vector<double> > > pooled_pattern_coherence;
    if (write_sites) {
        pooled_pattern_coherence.assign(
            analysed.size(),
            vector<vector<double> >(formula_specs.size(), vector<double>(aln->size(), 0.0)));
    }

    for (size_t category_index = 0; category_index < rate_categories.size(); category_index++) {
        const SatuTeRateCategory &category = rate_categories[category_index];

        satute::SatuTeDirectionalCache directional_cache;
        if (use_directional_cache) {
            satute::buildDirectionalCache(
                sweep_root, aln, model, nstates, category.rate, num_threads,
                node_bound, max_degree, directional_cache);
        }

        vector<double> category_pattern_coherence;
        for (size_t branch_index = 0; branch_index < analysed.size(); branch_index++) {
            const AnalysedBranch &entry = analysed[branch_index];
            double effective_length = entry.annotated_branch->length * category.rate;

            // The directed singular basis and the Frobenius information
            // fraction depend on the effective length, so both are rebuilt per
            // (branch, category) and shared by the two formulas.
            SatuTeDirectedBasis directed_basis;
            double directed_information = numeric_limits<double>::quiet_NaN();
            if (directed_nonreversible) {
                formulas::buildDirectedBasis(model, nstates, effective_length, directed_basis);
                directed_information = information::computeDirectedInformationFraction(
                    model, nstates, effective_length, entry.sibling_length * category.rate);
            }

            for (size_t formula_index = 0; formula_index < formula_specs.size(); formula_index++) {
                const SatuTeFormulaSpec &formula = formula_specs[formula_index];
                SatuTeFormulaSpec active_formula = formula;
                vector<double> weights;
                vector<double> left_projection;
                vector<double> right_projection;
                vector<double> mode_sigma;
                if (directed_nonreversible) {
                    active_formula.modes = formula.eigenvalue_weighted
                        ? formulas::allModeOffsets(directed_basis)
                        : formulas::dominantModeOffsets(directed_basis);
                    weights = formulas::buildSigmaWeights(
                        directed_basis, active_formula.modes, formula.eigenvalue_weighted);
                    for (size_t i = 0; i < active_formula.modes.size(); i++) {
                        int offset = active_formula.modes[i];
                        mode_sigma.push_back(directed_basis.sigma[offset]);
                        for (int state = 0; state < nstates; state++) {
                            left_projection.push_back(
                                directed_basis.left_projection[(size_t)offset * nstates + state]);
                            right_projection.push_back(
                                directed_basis.right_projection[(size_t)offset * nstates + state]);
                        }
                    }
                } else {
                    weights = formulas::buildFormulaWeights(eval, formula, effective_length);
                }

                SatuTeBranchResult base = results::makeBaseResult(
                    entry.branch_id,
                    entry.left_taxa,
                    entry.right_taxa,
                    entry.annotated_branch->length,
                    effective_length,
                    alpha,
                    formula.name,
                    category.label,
                    format::formatDouble(category.rate),
                    format::formatIntVector(active_formula.modes),
                    // Without a fixed eigenbasis the spectral column carries
                    // the canonical correlations sigma_k of this branch.
                    directed_nonreversible
                        ? format::formatVector(mode_sigma)
                        : format::eigenvalueString(eval, formula.modes),
                    format::formatVector(weights),
                    entry.label,
                    entry.split);
                information::setSaturationScale(
                    base,
                    directed_nonreversible
                        ? directed_information
                        : information::computeInformationFraction(eval, all_nonzero_modes, effective_length));

                SatuTeDirectedContext directed_context;
                if (directed_nonreversible) {
                    directed_context.left_projection = left_projection.data();
                    directed_context.right_projection = right_projection.data();
                    directed_context.parent_of = &parent_of;
                    directed_context.root_node_id = root->id;
                }

                category_stats[branch_index][formula_index].push_back(
                    computeCategoryStatistic(
                        entry.left_node,
                        entry.right_node,
                        entry.right_node,
                        entry.left_node,
                        aln,
                        model,
                        nstates,
                        eval,
                        active_formula,
                        weights,
                        pattern_cat_ptr,
                        category.pattern_cat,
                        category.rate,
                        num_threads,
                        enhanced_weights,
                        base,
                        directional_cache.ready ? &directional_cache : nullptr,
                        write_sites ? &category_pattern_coherence : nullptr,
                        directed_nonreversible ? &directed_context : nullptr));

                if (write_sites) {
                    vector<double> &target = pooled_pattern_coherence[branch_index][formula_index];
                    for (size_t ptn = 0; ptn < target.size() &&
                                         ptn < category_pattern_coherence.size(); ptn++)
                        target[ptn] += category_pattern_coherence[ptn];
                }
            }
        }
    }

    // Reduction and emission, branch-major exactly as before.
    for (size_t branch_index = 0; branch_index < analysed.size(); branch_index++) {
        const AnalysedBranch &entry = analysed[branch_index];
        const int branch_id = entry.branch_id;
        const int left_taxa = entry.left_taxa;
        const int right_taxa = entry.right_taxa;
        const string &label = entry.label;
        const string &split = entry.split;
        const double pooled_information_fraction = entry.pooled_information_fraction;
        Neighbor *annotated_branch = entry.annotated_branch;

        for (size_t formula_index = 0; formula_index < formula_specs.size(); formula_index++) {
            const SatuTeFormulaSpec &formula = formula_specs[formula_index];
            const vector<SatuTeBranchResult> &category_results =
                category_stats[branch_index][formula_index];

            SatuTeBranchResult pooled_base = results::makeBaseResult(
                branch_id,
                left_taxa,
                right_taxa,
                annotated_branch->length,
                annotated_branch->length,
                alpha,
                formula.name,
                "pooled",
                "NA",
                // Directed mode sets and canonical correlations are
                // category-specific, so pooled rows carry no single set.
                directed_nonreversible ? "NA" : format::formatIntVector(formula.modes),
                directed_nonreversible ? "NA" : format::eigenvalueString(eval, formula.modes),
                "NA",
                label,
                split);
            information::setSaturationScale(pooled_base, pooled_information_fraction);

            SatuTeBranchResult pooled = (rate_categories.size() == 1 && rate_categories[0].pattern_cat < 0)
                ? category_results[0]
                : poolCategoryResults(category_results, pooled_base);
            pooled.rate_category = "pooled";
            pooled.rate_multiplier = "NA";
            pooled.effective_length = pooled.length;
            information::setSaturationScale(pooled, pooled_information_fraction);
            results.push_back(pooled);

            if (write_sites) {
                SatuTeSiteSeries series;
                series.branch_id = branch_id;
                series.formula = formula.name;
                series.split = split;
                series.coherence = pooled.coherence;
                series.p_value = pooled.p_value;
                series.valid_sites = pooled.valid_sites;
                series.pattern_coherence =
                    pooled_pattern_coherence[branch_index][formula_index];
                site_series.push_back(series);
            }

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
            annotated_branch->putAttr("satDominantFDR", format::formatDouble(result.fdr_by));
            annotated_branch->putAttr("satDominantFDRDecision", result.decision_fdr);
            continue;
        }
        if (result.formula != "eigenvalue_weighted")
            continue;

        annotated_branch->putAttr("satFormula", result.formula);
        annotated_branch->putAttr("satC", format::formatDouble(result.coherence));
        annotated_branch->putAttr("satSE", format::formatDouble(result.se));
        annotated_branch->putAttr("satZ", format::formatDouble(result.z_score));
        annotated_branch->putAttr("satP", format::formatDouble(result.p_value));
        annotated_branch->putAttr("satFDR", format::formatDouble(result.fdr_by));
        annotated_branch->putAttr("satFDRDecision", result.decision_fdr);
        annotated_branch->putAttr("satInfo", format::formatDouble(result.information_fraction));
        annotated_branch->putAttr("satIndex", format::formatDouble(result.saturation_index));
        annotated_branch->putAttr("sat", result.decision);
        annotated_branch->putAttr("satBonf", result.decision_bonferroni);
    }

    ReportWriter().write(*this, prefix, results);
    if (write_sites) {
        // Keep the branches closest to saturation. The whole-tree export is
        // impractical on large trees (~507 MB for the 1871-taxon Tree of
        // Life), and the branches worth windowing are the marginal ones.
        if (sites_max > 0 && (int)site_series.size() > sites_max) {
            stable_sort(
                site_series.begin(),
                site_series.end(),
                [](const SatuTeSiteSeries &a, const SatuTeSiteSeries &b) {
                    return a.p_value > b.p_value;
                });
            site_series.resize(sites_max);
            cout << "SatuTe per-site export limited to " << sites_max
                 << " branches with the highest satP" << endl;
        }
        ReportWriter().writeSiteSeries(*this, prefix, site_series);
    }
}
