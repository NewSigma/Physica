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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.cuh"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<>;

namespace {
    void rvalue_host_device() {
        constexpr int N = 128;
        using T = float32;
        const auto A = MatrixND<T>::random_uniform<RandomSource>(N, N);
        const auto x = VectorND<T>::random_uniform<RandomSource>(N);
        VectorND<T> y = A * exp(x);
        auto dy = device_obj<VectorND<T>>(A.toDevice() * exp(x.toDevice())).toHost();
        expect<RandomSource>(vectorNear(y, dy, 1E-6));

        y = A.transpose() * exp(x);
        dy = device_obj<VectorND<T>>(A.toDevice().transpose() * exp(x.toDevice())).toHost();
        expect<RandomSource>(vectorNear(y, dy, 1E-6));
    }

    void cublas() {
        constexpr int N = 128;
        using T = float32;
        const auto A = DenseMatrix<T, MatrixMajor::Row>::random_uniform<RandomSource>(N, N);
        const auto x = VectorND<T>::random_uniform<RandomSource>(N);
        const VectorND<T> answer = A * x;
        device_obj<VectorND<T>> dy(N);
        (A.toDevice() * x.toDevice()).assign_cublas(dy);
        auto result = dy.toHost();
        expect<RandomSource>(vectorNear(result, answer, 4UL));
    }
}

int main() {
    rvalue_host_device();
    cublas();
    return 0;
}
