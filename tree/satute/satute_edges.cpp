/***************************************************************************
 *   Copyright (C) 2026 by                                                *
 *   IQ-TREE developers                                                    *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 ***************************************************************************/

#include "satute_edges.h"

#include "../phylotree.h"

#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>

using namespace std;

namespace satute {
namespace edges {

set<int> readEdgeFile(const char *filename) {
    set<int> edges;
    if (filename == nullptr || strlen(filename) == 0)
        return edges;

    ifstream in(filename);
    if (!in.is_open())
        outError("Cannot read SatuTe branch ID file ", filename);

    string line;
    while (getline(in, line)) {
        size_t comment = line.find('#');
        if (comment != string::npos)
            line.erase(comment);

        stringstream ss(line);
        string token;
        while (ss >> token) {
            char *end = nullptr;
            long id = strtol(token.c_str(), &end, 10);
            if (end == token.c_str() || *end != '\0')
                outError("Invalid SatuTe branch ID in file ", filename);
            if (id < numeric_limits<int>::min() || id > numeric_limits<int>::max())
                outError("SatuTe branch ID out of integer range in file ", filename);
            edges.insert(id);
        }
    }

    if (edges.empty())
        outError("SatuTe branch ID file is empty");

    return edges;
}

} // namespace edges
} // namespace satute
