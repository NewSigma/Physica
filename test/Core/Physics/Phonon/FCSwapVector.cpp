/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Physics/Phonon/PhononSolverImpl/FCSwapVector.h"
#include "Test.h"

using namespace Physica;
using ScalarType = float64;
using VectorType = FCSwapVector<ScalarType>;

namespace {
    void testIndex1D5D(size_t numDOF, const Index3D& superSize, size_t index1D) {
        const auto index5D = VectorType::index1DTo5D(numDOF, superSize, index1D);
        const auto index1D_1 = VectorType::index5DTo1D(numDOF, superSize, index5D);
        expect(index1D == index1D_1);
    }
}

int main() {
    const size_t numDOF = 12;
    const Index3D superSize{3, 4, 5};
    testIndex1D5D(numDOF, superSize, 173);
    testIndex1D5D(numDOF, superSize, 997);
    testIndex1D5D(36, {8, 8, 1}, 1296);
    return 0;
}
