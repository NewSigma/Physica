/*
 * Copyright 2022-2025 Weibo He.
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
#include "Physica/Core/ML/Clustering/AP.h"
#include "Test.h"

using namespace Physica;

namespace {
    bool testCluster(const std::forward_list<size_t>& cluster) {
        const size_t first = cluster.front();
        for (auto i : cluster) {
            if (first < 3)
                expect(i < 3);
            else
                expect(i >= 3);
        }
        return true;
    }
}

int main() {
    DenseMatrix<float64> points{{5, 0}, {5, 1}, {5, -1}, {-5, 0}, {-5, 1}, {-5, -1}};
    DenseSymmMatrix<float64> similar(points.getCol());
    for (size_t i = 0; i < similar.getOrder(); ++i) {
        similar[i, i] = float64(-1);
        for (size_t j = i + 1; j < similar.getOrder(); ++j) {
            similar[i, j] = -(points.col(i) - points.col(j)).norm();
        }
    }

    AP<float64> ap(similar, 0.5, 10, 100);
    auto exemplars = ap.getExemplars();
    expect(exemplars.size() == 2);
    {
        auto ite = exemplars.begin();
        testCluster(ap.getCluster(*ite));

        ++ite;
        testCluster(ap.getCluster(*ite));
    }
    return 0;
}
