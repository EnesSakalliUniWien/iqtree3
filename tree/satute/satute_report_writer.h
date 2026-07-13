/***************************************************************************
 *   Copyright (C) 2026 by                                                *
 *   IQ-TREE developers                                                    *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 ***************************************************************************/

#ifndef SATUTE_REPORT_WRITER_H
#define SATUTE_REPORT_WRITER_H

#include "satute_types.h"
#include "../phylotree.h"

#include <string>
#include <vector>

namespace satute {

class ReportWriter {
public:
    void write(
        PhyloTree &tree,
        const std::string &prefix,
        const std::vector<SatuTeBranchResult> &results) const;
};

} // namespace satute

#endif
