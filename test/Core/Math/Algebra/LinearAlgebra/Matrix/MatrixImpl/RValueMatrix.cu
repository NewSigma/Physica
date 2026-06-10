/*
 * Copyright 2025-2026 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.cuh"
#include "Test.h"

using namespace Physica;
using MatrixType = MatrixND<float32>;

namespace {
    void sum() {
        const MatrixType A = MatrixType::random_uniform<Random<>>(16, 16);
        const auto d_A = A.toDeviceAsync();
        expect(scalarNear(d_A.diag().sum(), A.diag().sum(), 1E-6));
    }

    void minorDiag() {
        const MatrixType A = MatrixType::random_uniform<Random<>>(4, 4);
        const auto d_A = A.toDeviceAsync();
        int shift = Array<int>{1, 2, 3, -1, -2, -3, 0}.template select<Random<>>();

        device_obj<VectorND<float32>> result = (-d_A).diag(shift);
        expect((-A).diag(shift) == result.toHost());
    }
}

static_assert(std::same_as<decltype(device_obj<MatrixND<float32>>{}.transpose().transpose()), device_obj<MatrixND<float32>>&&>, "[Error]: transpose-transpose should yield it self");

int main() {
    sum();
    minorDiag();
    return 0;
}
