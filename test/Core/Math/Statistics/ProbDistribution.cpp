/*
 * Copyright 2026 Weibo He.
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
#include "Physica/Core/Math/Statistics/ProbDistribution.h"
#include "Physica/Core/Math/Statistics/ProbDistribution2D.h"
#include "Test.h"

using namespace Physica;
using T = float64;
using RandomSource = Random<MCG>;

int main() {
    // Test that getMin(), getMax(), and calcNumSample() are helpful diagnostics when data lies outside the expected range
    const auto data = VectorND<T>::random_uniform<RandomSource>(8);
    {
        ProbDistribution<T> dist(2, 3, 10);
        dist.sample(data);
        expect(dist.getMinimum() == data.min());
        expect(dist.getMaximum() == data.max());
        expect(dist.calcNumSample() == 0);
    }
    {
        ProbDistribution2D<T> dist(2, 3, 2, 3, 8, 8);
        for (auto elem : data)
            dist.sample(elem, elem);
        expect(dist.getMinX() == data.min());
        expect(dist.getMaxX() == data.max());
        expect(dist.getMinY() == data.min());
        expect(dist.getMaxY() == data.max());
        expect(dist.calcNumSample() == 0);
    }
    return 0;
}
