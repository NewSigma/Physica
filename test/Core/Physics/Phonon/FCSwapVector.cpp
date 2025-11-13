/*
 * Copyright 2024-2025 Weibo He.
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
#include "Physica/Core/Physics/Phonon/PhononSolverImpl/FCSwapVector.h"

using namespace Physica;
using ScalarType = float64;
using VectorType = FCSwapVector<ScalarType>;

namespace {
    bool testIndex1D5D(size_t numDOF, const Index3D& superSize, size_t index1D) {
        const auto index5D = VectorType::index1DTo5D(numDOF, superSize, index1D);
        const auto index1D_1 = VectorType::index5DTo1D(numDOF, superSize, index5D);
        return index1D == index1D_1;
    }
}

int main() {
    const size_t numDOF = 12;
    const Index3D superSize{3, 4, 5};
    if (!testIndex1D5D(numDOF, superSize, 173))
        return 1;
    if (!testIndex1D5D(numDOF, superSize, 997))
        return 1;
    if (!testIndex1D5D(36, {8, 8, 1}, 1296))
        return 1;
    return 0;
}
