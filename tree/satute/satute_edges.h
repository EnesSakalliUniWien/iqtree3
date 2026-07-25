/***************************************************************************
 *   Copyright (C) 2026 by                                                *
 *   IQ-TREE developers                                                    *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 ***************************************************************************/

#ifndef SATUTE_EDGES_H
#define SATUTE_EDGES_H

#include <set>

namespace satute {
namespace edges {

/** Parse the optional --satute-edges file restricting analysis to branch IDs. */
std::set<int> readEdgeFile(const char *filename);

} // namespace edges
} // namespace satute

#endif
