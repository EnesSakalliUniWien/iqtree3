/***************************************************************************
 *   Copyright (C) 2026 by                                                *
 *   IQ-TREE developers                                                    *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 ***************************************************************************/

#ifndef SATUTE_TYPES_H
#define SATUTE_TYPES_H

#include <string>
#include <vector>

namespace satute {

struct SatuTeFormulaSpec {
    std::string name;
    std::vector<int> modes;
    bool eigenvalue_weighted;
};

struct SatuTeRateCategory {
    std::string label;
    int pattern_cat;
    double rate;
    double proportion;
};

struct SatuTeBranchResult {
    int id;
    int left_taxa;
    int right_taxa;
    int valid_sites;
    int skipped_sites;
    double length;
    double effective_length;
    double information_fraction;
    double saturation_index;
    double coherence;
    double variance;
    double se;
    double z_score;
    double p_value;
    double alpha;
    double alpha_adjusted;
    std::string formula;
    std::string rate_category;
    std::string rate_multiplier;
    int rate_category_sites;
    std::string modes;
    std::string eigenvalues;
    std::string weights;
    std::string decision;
    std::string decision_bonferroni;
    std::string label;
    std::string split;
};


} // namespace satute

#endif
