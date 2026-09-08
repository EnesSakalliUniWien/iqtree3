/***************************************************************************
 *   Copyright (C) 2026 by                                                *
 *   IQ-TREE developers                                                    *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 ***************************************************************************/

#ifndef SATUTE_TOPOLOGY_H
#define SATUTE_TOPOLOGY_H

#include "../phylotree.h"

#include <string>
#include <vector>

namespace satute {
namespace topology {

/** Tip counts and split labels for a branch, used to describe each result row.
 * The __root__ leaf of a rooted tree is an attachment point, not a taxon, and
 * is excluded from both. */
int countTipsAwayFrom(Node *node, Node *dad);
std::string splitLabel(Node *left, Node *right);

/**
 * Parent relation of the rooted tree: parent_of[node id] = parent node id,
 * -1 for the root itself. Needed by the direction-aware non-reversible
 * recursion to decide which end of every edge is the true parent.
 */
void fillParentMap(Node *node, Node *dad, std::vector<int> &parent_of);

} // namespace topology
} // namespace satute

#endif
