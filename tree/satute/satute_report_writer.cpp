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
#include "satute_format.h"

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
    if (!out.is_open())
        outError("Cannot write SatuTe statistics file ", filename);
    out << "# SatuTe branch saturation statistics" << endl
        << "# This file can be read in MS Excel or in R with command:" << endl
        << "#   tab=read.table('" << filename << "',header=TRUE)" << endl
        << "# Columns are tab-separated with following meaning:" << endl
        << "#   ID: Branch ID" << endl
        << "#   Formula: dominant or eigenvalue_weighted" << endl
        << "#   RateCategory: pooled for the test, or an IQ-TREE maximum-posterior category diagnostic" << endl
        << "#   RateMultiplier: IQ-TREE category rate used to rescale branch lengths; NA for pooled rows" << endl
        << "#   RateSites: Number of sites assigned by IQ-TREE to this rate category or pooled row" << endl
        << "#   LeftTaxa, RightTaxa: Number of taxa on each side of the branch" << endl
        << "#   ValidSites: Alignment sites used after skipping invalid likelihood ratios" << endl
        << "#   SkippedSites: Alignment sites skipped because likelihood ratios were invalid" << endl
        << "#   satC: Mean SatuTe branch coherence statistic" << endl
        << "#   satVar: Variance estimate used to studentize satC" << endl
        << "#     homogeneous tests use factorized marginal covariances; pooled rate mixtures use the sample variance of the centred integrated site scores" << endl
        << "#     hard-category moment estimates are diagnostic" << endl
        << "#   satSE: Standard error of pooled satC; diagnostic scale for hard categories" << endl
        << "#   satZ: Pooled one-sided normal test statistic; diagnostic standardized coherence for hard categories" << endl
        << "#   satP: One-sided p-value; NA for hard-category diagnostics, which do not define independent null samples" << endl
        << "#   Alpha: Unadjusted significance level" << endl
        << "#   AlphaTaxonBonf: Alpha adjusted by LeftTaxa*RightTaxa" << endl
        << "#   Decision: informative if satP <= Alpha, otherwise saturated; diagnostic for hard categories" << endl
        << "#   DecisionTaxonBonf: informative if satP <= AlphaTaxonBonf" << endl
        << "#   FDR_BY: Benjamini-Yekutieli adjusted p-value for pooled rows; dominant and eigenvalue_weighted are adjusted as separate branch families" << endl
        << "#   DecisionFDR: informative if FDR_BY <= Alpha; not_tested for rate-category rows" << endl
        << "#   Label: Existing branch label" << endl
        << "#   Length: Original branch length" << endl
        << "#   EffectiveLength: Branch length after category-rate rescaling" << endl
        << "#   InformationFraction: Model-based non-stationary spectral energy I(t)/I(0), integrated over rate proportions for pooled rows" << endl
        << "#   SaturationIndex: Monotonic model-based scale 1-InformationFraction; 0 means no spectral decay and 1 means complete model saturation" << endl
        << "#   Split: Comma-separated taxa on the smaller side of the branch" << endl
        << "#   Modes, Eigenvalues, Weights: spectral modes used by the formula; pooled mixture rows use NA weights" << endl
        << "# The annotated .sat.tree and .sat.tree.nex files use the pooled eigenvalue_weighted row and record both formula-specific FDR results" << endl
        << "ID\tFormula\tRateCategory\tRateMultiplier\tRateSites\tLeftTaxa\tRightTaxa\tValidSites\tSkippedSites\tsatC\tsatVar\tsatSE\tsatZ\tsatP\tAlpha\tAlphaTaxonBonf\tDecision\tDecisionTaxonBonf\tFDR_BY\tDecisionFDR\tLabel\tLength\tEffectiveLength\tInformationFraction\tSaturationIndex\tSplit\tModes\tEigenvalues\tWeights" << endl;

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
            << '\t' << format::formatDouble(result.coherence)
            << '\t' << format::formatDouble(result.variance)
            << '\t' << format::formatDouble(result.se)
            << '\t' << format::formatDouble(result.z_score)
            << '\t' << format::formatDouble(result.p_value)
            << '\t' << format::formatDouble(result.alpha)
            << '\t' << format::formatDouble(result.alpha_adjusted)
            << '\t' << result.decision
            << '\t' << result.decision_bonferroni
            << '\t' << format::formatDouble(result.fdr_by)
            << '\t' << result.decision_fdr
            << '\t' << result.label
            << '\t' << format::formatDouble(result.length)
            << '\t' << format::formatDouble(result.effective_length)
            << '\t' << format::formatDouble(result.information_fraction)
            << '\t' << format::formatDouble(result.saturation_index)
            << '\t' << result.split
            << '\t' << result.modes
            << '\t' << result.eigenvalues
            << '\t' << result.weights
            << endl;
    }
    out.close();
    cout << "SatuTe statistics per branch printed to " << filename << endl;
}

void ReportWriter::writeSiteSeries(
    PhyloTree &tree,
    const string &prefix,
    const vector<SatuTeSiteSeries> &series) const
{
    if (series.empty())
        return;

    string filename = prefix + ".sat.site";
    ofstream out(filename.c_str());
    if (!out.is_open())
        outError("Cannot write SatuTe per-site file ", filename);

    Alignment *aln = tree.aln;
    size_t nsite = aln->getNSite();

    out << "# SatuTe per-site contributions to satC" << endl
        << "# One row per branch, formula and alignment site." << endl
        << "# Columns are tab-separated with following meaning:" << endl
        << "#   ID: Branch ID (matches the ID column in .sat.stat)" << endl
        << "#   Formula: dominant or eigenvalue_weighted" << endl
        << "#   Site: 1-based alignment column" << endl
        << "#   SiteCoherence: this site's centred contribution, homogeneous or integrated over rate categories;" << endl
        << "#     the frequency-weighted mean of these over valid sites equals" << endl
        << "#     the branch's pooled satC in .sat.stat" << endl
        << "# NOTE: sites sharing an alignment pattern carry identical values." << endl
        << "# Branch splits are not repeated per row -- look up the ID in .sat.stat." << endl
        << "ID\tFormula\tSite\tSiteCoherence" << endl;

    for (size_t i = 0; i < series.size(); i++) {
        const SatuTeSiteSeries &entry = series[i];
        for (size_t site = 0; site < nsite; site++) {
            int ptn = aln->getPatternID(site);
            if (ptn < 0 || (size_t)ptn >= entry.pattern_coherence.size())
                continue;
            out << entry.branch_id
                << '\t' << entry.formula
                << '\t' << (site + 1)
                << '\t' << format::formatDouble(entry.pattern_coherence[ptn])
                << endl;
        }
    }
    out.close();
    cout << "SatuTe per-site contributions printed to " << filename << endl;
}

} // namespace satute
