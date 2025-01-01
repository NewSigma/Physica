/*
 * Copyright 2022 Weibo He.
 *
 * This file is part of Physica.
 *
 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
 */
#include <iostream>
#include "Physica/Core/AI/Clustering/AP.h"

using namespace Physica;
using ScalarType = float64;

bool testCluster(const std::forward_list<size_t>& cluster) {
    const size_t first = cluster.front();
    for (auto ite = cluster.cbegin(); ite != cluster.cend(); ++ite) {
        if (first < 3 && *ite >= 3)
            return false;
        if (first >= 3 && *ite < 3)
            return false;
    }
    return true;
}

int main() {
    DenseMatrix<ScalarType> points{{5, 0}, {5, 1}, {5, -1}, {-5, 0}, {-5, 1}, {-5, -1}};
    DenseSymmMatrix<ScalarType> similar(points.getCol());
    for (size_t i = 0; i < similar.getOrder(); ++i) {
        similar(i, i) = ScalarType(-1);
        for (size_t j = i + 1; j < similar.getOrder(); ++j) {
            similar(i, j) = -(points.col(i) - points.col(j)).norm();
        }
    }

    AP<ScalarType> ap(std::move(similar), 0.5, 10, 100);
    auto exemplars = ap.getExemplars();
    if (exemplars.size() != 2)
        return 1;
    {
        auto ite = exemplars.begin();
        if (!testCluster(ap.getCluster(*ite)))
            return 1;
        ++ite;
        if (!testCluster(ap.getCluster(*ite)))
            return 1;
    }
    return 0;
}
