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

namespace satute {
namespace topology {

/** Tip counts and split labels for a branch, used to describe each result row. */
int countTipsAwayFrom(Node *node, Node *dad);
std::string splitLabel(Node *left, Node *right);

} // namespace topology
} // namespace satute

#endif
