/***************************************************************************
 *   Copyright (C) 2026 by                                                *
 *   IQ-TREE developers                                                    *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 ***************************************************************************/

#include "satute_format.h"

#include <cmath>
#include <iomanip>
#include <sstream>

using namespace std;

namespace satute {
namespace format {

string formatDouble(double value) {
    if (!isfinite(value))
        return "NA";
    stringstream ss;
    ss << setprecision(10) << value;
    return ss.str();
}

string formatVector(const vector<double> &values) {
    stringstream ss;
    ss << setprecision(10);
    for (size_t i = 0; i < values.size(); i++) {
        if (i > 0)
            ss << ',';
        ss << values[i];
    }
    return ss.str();
}

string formatIntVector(const vector<int> &values) {
    stringstream ss;
    for (size_t i = 0; i < values.size(); i++) {
        if (i > 0)
            ss << ',';
        ss << values[i];
    }
    return ss.str();
}

string eigenvalueString(const double *eval, const vector<int> &modes) {
    vector<double> values;
    for (size_t i = 0; i < modes.size(); i++)
        values.push_back(eval[modes[i]]);
    return formatVector(values);
}

} // namespace format
} // namespace satute
