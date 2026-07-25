/***************************************************************************
 *   Copyright (C) 2026 by                                                *
 *   IQ-TREE developers                                                    *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 ***************************************************************************/

#ifndef SATUTE_FORMAT_H
#define SATUTE_FORMAT_H

#include <string>
#include <vector>

namespace satute {
namespace format {

/** Numeric-to-text helpers shared by the report writer and the driver. */
std::string formatDouble(double value);
std::string formatVector(const std::vector<double> &values);
std::string formatIntVector(const std::vector<int> &values);
std::string eigenvalueString(const double *eval, const std::vector<int> &modes);

} // namespace format
} // namespace satute

#endif
