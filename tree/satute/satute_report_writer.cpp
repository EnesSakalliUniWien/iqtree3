/***************************************************************************
 *   Copyright (C) 2026 by                                                *
 *   IQ-TREE developers                                                    *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 ***************************************************************************/

#include "satute_report_writer.h"
#include "satute_support.h"

#include <fstream>
#include <iostream>

using namespace std;

namespace satute {

void ReportWriter::write(
    PhyloTree &tree,
    const string &prefix,
    const vector<SatuTeBranchResult> &results) const
{
    string out_prefix = prefix;
    string filename = out_prefix + ".sat.stat";
    string str = out_prefix + ".sat.tree";
    tree.printTree(str.c_str(), WT_BR_LEN + WT_NEWLINE);
    cout << "Tree with SatuTe annotations written to " << str << endl;

    str = out_prefix + ".sat.tree.nex";
    tree.printNexus(str, WT_BR_LEN, "See " + filename + " for SatuTe branch saturation statistics. "
               "This file is best viewed in FigTree.");
    cout << "Annotated tree (best viewed in FigTree) written to " << str << endl;

    str = out_prefix + ".sat.branch";
    tree.printTree(str.c_str(), WT_BR_LEN + WT_INT_NODE + WT_NEWLINE);
    cout << "Tree with SatuTe branch IDs written to " << str << endl;

    ofstream out(filename.c_str());
    out << "# SatuTe branch saturation statistics" << endl
        << "# This file can be read in MS Excel or in R with command:" << endl
        << "#   tab=read.table('" << filename << "',header=TRUE)" << endl
        << "# Columns are tab-separated with following meaning:" << endl
        << "#   ID: Branch ID" << endl
        << "#   Formula: dominant, eigenvalue_weighted, or mixture_likelihood_weighted" << endl
        << "#   mixture_likelihood_weighted: soft null-category integration with one common log-weight shift; zero-rate invariant categories enter only the infinite-branch null denominator" << endl
        << "#   RateCategory: pooled for all sites, or the IQ-TREE maximum-posterior pattern-rate category" << endl
        << "#   RateMultiplier: IQ-TREE category rate used to rescale branch lengths; NA for pooled rows" << endl
        << "#   RateSites: Number of sites assigned by IQ-TREE to this rate category or pooled row" << endl
        << "#   LeftTaxa, RightTaxa: Number of taxa on each side of the branch" << endl
        << "#   ValidSites: Alignment sites used after skipping invalid likelihood ratios" << endl
        << "#   SkippedSites: Alignment sites skipped because likelihood ratios were invalid" << endl
        << "#   satC: Mean SatuTe branch coherence statistic" << endl
        << "#   satVar: Variance estimate used to studentize satC" << endl
        << "#     dominant/eigenvalue_weighted use the factorized saturated-null variance from left/right marginal second moments" << endl
        << "#     mixture_likelihood_weighted uses the ordinary sample variance of its final soft-mixture site scores" << endl
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
        << "#   InformationFraction: Model-based non-stationary spectral energy I(t)/I(0), integrated over rate proportions for pooled rows" << endl
        << "#   SaturationIndex: Monotonic model-based scale 1-InformationFraction; 0 means no spectral decay and 1 means complete model saturation" << endl
        << "#   Split: Comma-separated taxa on the smaller side of the branch" << endl
        << "#   Modes, Eigenvalues, Weights: spectral modes used by the formula; pooled hard-category rows use NA weights, while the soft-mixture row reports its common global log-weight shift" << endl
        << "# The annotated .sat.tree and .sat.tree.nex files use the pooled eigenvalue_weighted row and record satFormula=eigenvalue_weighted" << endl
        << "ID\tFormula\tRateCategory\tRateMultiplier\tRateSites\tLeftTaxa\tRightTaxa\tValidSites\tSkippedSites\tsatC\tsatVar\tsatSE\tsatZ\tsatP\tAlpha\tAlphaTaxonBonf\tDecision\tDecisionTaxonBonf\tLabel\tLength\tEffectiveLength\tInformationFraction\tSaturationIndex\tSplit\tModes\tEigenvalues\tWeights" << endl;

    for (size_t i = 0; i < results.size(); i++) {
        const SatuTeBranchResult &result = results[i];
        out << result.id
            << '\t' << result.formula
            << '\t' << result.rate_category
            << '\t' << result.rate_multiplier
            << '\t' << result.rate_category_sites
            << '\t' << result.left_taxa
            << '\t' << result.right_taxa
            << '\t' << result.valid_sites
            << '\t' << result.skipped_sites
            << '\t' << formatDouble(result.coherence)
            << '\t' << formatDouble(result.variance)
            << '\t' << formatDouble(result.se)
            << '\t' << formatDouble(result.z_score)
            << '\t' << formatDouble(result.p_value)
            << '\t' << formatDouble(result.alpha)
            << '\t' << formatDouble(result.alpha_adjusted)
            << '\t' << result.decision
            << '\t' << result.decision_bonferroni
            << '\t' << result.label
            << '\t' << formatDouble(result.length)
            << '\t' << formatDouble(result.effective_length)
            << '\t' << formatDouble(result.information_fraction)
            << '\t' << formatDouble(result.saturation_index)
            << '\t' << result.split
            << '\t' << result.modes
            << '\t' << result.eigenvalues
            << '\t' << result.weights
            << endl;
    }
    out.close();
    cout << "SatuTe statistics per branch printed to " << filename << endl;
}

} // namespace satute
