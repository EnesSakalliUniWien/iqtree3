/***************************************************************************
 *   Copyright (C) 2026 by                                                *
 *   IQ-TREE developers                                                    *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 ***************************************************************************/

#include "satute_topology.h"

#include <algorithm>
#include <cstring>
#include <sstream>
#include <vector>

using namespace std;

namespace satute {
namespace topology {

/** The __root__ leaf of a rooted tree is an attachment point, not a taxon.
 * ROOT_NAME never occurs as a real taxon name, so this is inert on the
 * unrooted reversible path. */
static bool isRootAttachment(Node *node) {
    return node->isLeaf() && node->name == ROOT_NAME;
}

int countTipsAwayFrom(Node *node, Node *dad) {
    if (node->isLeaf())
        return isRootAttachment(node) ? 0 : 1;

    int count = 0;
    FOR_NEIGHBOR_IT(node, dad, it) {
        count += countTipsAwayFrom((*it)->node, node);
    }
    return count;
}

void collectTipsAwayFrom(Node *node, Node *dad, vector<string> &tips) {
    if (node->isLeaf()) {
        if (!isRootAttachment(node))
            tips.push_back(node->name);
        return;
    }

    FOR_NEIGHBOR_IT(node, dad, it) {
        collectTipsAwayFrom((*it)->node, node, tips);
    }
}

void fillParentMap(Node *node, Node *dad, vector<int> &parent_of) {
    if (dad == nullptr)
        parent_of[node->id] = -1;

    FOR_NEIGHBOR_IT(node, dad, it) {
        parent_of[(*it)->node->id] = node->id;
        fillParentMap((*it)->node, node, parent_of);
    }
}

string splitLabel(Node *left, Node *right) {
    vector<string> left_tips;
    vector<string> right_tips;
    collectTipsAwayFrom(left, right, left_tips);
    collectTipsAwayFrom(right, left, right_tips);

    vector<string> &chosen = (left_tips.size() <= right_tips.size()) ? left_tips : right_tips;
    sort(chosen.begin(), chosen.end());

    stringstream ss;
    for (size_t i = 0; i < chosen.size(); i++) {
        if (i > 0)
            ss << ',';
        ss << chosen[i];
    }
    return ss.str();
}

} // namespace topology
} // namespace satute
