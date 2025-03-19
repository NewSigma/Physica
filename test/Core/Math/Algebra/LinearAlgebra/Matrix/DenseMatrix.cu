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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.cuh"

using namespace Physica;
using RandomType = Random<MT19937, 10000>;

int main() {
    {
        using MatrixType = DenseMatrix<float32>;
        const MatrixType A = MatrixType::random_uniform<RandomType>(16, 16);
        const auto d_A = A.toDevice();
        const MatrixType B = d_A.toHost();
        if (A.asArray() != B.asArray())
            return 1;
    }
    {
        using MatrixType = DenseMatrix<float32, MatrixOption::Col | MatrixOption::Element>;
        using DeviceMatrix = MatrixType::device_obj_type;
        const MatrixType A = MatrixType::random_uniform<RandomType>(3, 4);
        const MatrixType B = A.transpose();
        const MatrixType answer = A * B;

        const auto d_A = A.toDevice();
        const auto d_B = B.toDevice();
        const DeviceMatrix d_result = d_A * d_B;
        CUDAContext::getInstance().wait();
        const auto result = d_result.toHost();
        if (!matrixNear(answer, result, 1E-6))
            return 1;
    }
    return 0;
}
