/***************************************************************************
 *   Copyright (C) 2026 by                                                *
 *   IQ-TREE developers                                                    *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 ***************************************************************************/

#ifndef SATUTE_SUPPORT_H
#define SATUTE_SUPPORT_H

#include "satute_types.h"
#include "../phylotree.h"

#include <set>
#include <string>
#include <vector>

namespace satute {

std::string formatDouble(double value);
std::string formatVector(const std::vector<double> &values);
std::string formatIntVector(const std::vector<int> &values);
int countTipsAwayFrom(Node *node, Node *dad);
std::string splitLabel(Node *left, Node *right);
std::set<int> readEdgeFile(const char *filename);
int findZeroEigenIndex(const double *eval, int nstates);
std::vector<int> findNonZeroModes(const double *eval, int nstates, int zero_index);
std::vector<SatuTeFormulaSpec> buildFormulaSpecs(const double *eval, int nstates, int zero_index);
std::vector<double> buildFormulaWeights(
    const double *eval,
    const SatuTeFormulaSpec &formula,
    double effective_branch_length);
std::string eigenvalueString(const double *eval, const std::vector<int> &modes);
double computeInformationFraction(
    const double *eval,
    const std::vector<int> &modes,
    double effective_length);
double computeMixtureInformationFraction(
    const double *eval,
    const std::vector<int> &modes,
    const std::vector<SatuTeRateCategory> &categories,
    double branch_length);
double computeMixtureLogWeightShift(
    const double *eval,
    const std::vector<int> &modes,
    const std::vector<SatuTeRateCategory> &categories,
    double branch_length);
void setSaturationScale(SatuTeBranchResult &result, double information_fraction);
SatuTeBranchResult makeBaseResult(
    int branch_id,
    int left_taxa,
    int right_taxa,
    double branch_length,
    double effective_length,
    double alpha,
    const std::string &formula,
    const std::string &rate_category,
    const std::string &rate_multiplier,
    const std::string &modes,
    const std::string &eigenvalues,
    const std::string &weights,
    const std::string &label,
    const std::string &split);
std::vector<SatuTeRateCategory> buildRateCategories(
    RateHeterogeneity *site_rate,
    DoubleVector &pattern_rates,
    IntVector &pattern_cat);

} // namespace satute

#endif
