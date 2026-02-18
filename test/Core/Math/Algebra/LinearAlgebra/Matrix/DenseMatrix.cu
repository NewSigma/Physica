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
#include "Test.h"

using namespace Physica;
using RandomSource = Random<MT19937, 10000>;

namespace {
    void hostDeviceCopy() {
        using MatrixType = DenseMatrix<float32>;
        const MatrixType A = MatrixType::random_uniform<RandomSource>(16, 16);
        const auto d_A = A.toDevice();
        const MatrixType B = d_A.toHost();
        expect(A.asArray() == B.asArray());
    }

    void deviceExprEval() {
        using MatrixType = DenseMatrix<float32, MatrixMajor::Col>;
        using DeviceMatrix = MatrixType::device_obj_type;
        const MatrixType A = MatrixType::random_uniform<RandomSource>(3, 4);
        const MatrixType B = A.transpose();
        const MatrixType answer = A * B;

        const auto d_A = A.toDevice();
        const auto d_B = B.toDevice();
        const DeviceMatrix d_result = d_A * d_B;
        CUDAContext::getInstance().wait();
        const auto result = d_result.toHost();
        expect(matrixNear(answer, result, 1E-6));
    }
    /**
     * A compact matrix is compact either in row or in column.
     */
    void CompactRowCol() noexcept {
        constexpr int Size = 8;
        using T = float32;
        auto& rng = RandomSource::getInstance();
        int r = std::uniform_int_distribution<>(0, Size - 1)(rng);
        int c = std::uniform_int_distribution<>(0, Size - 1)(rng);
        {
            using MatrixType = device_obj<DenseMatrix<T, MatrixMajor::Col>>;
            const auto x = MatrixType(Size, Size);
            expect(x.data_ptr(r, c) == x.col(c).data() + r);
        }
        {
            using MatrixType = device_obj<DenseMatrix<T, MatrixMajor::Row>>;
            const auto x = MatrixType(Size, Size);
            expect(x.data_ptr(r, c) == x.row(r).data() + c);
        }
    }
}

int main() {
    hostDeviceCopy();
    deviceExprEval();
    CompactRowCol();
    return 0;
}
