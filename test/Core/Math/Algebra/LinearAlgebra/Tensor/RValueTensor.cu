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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.cuh"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.cuh"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/DenseTensor.cuh"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Parallel/Executor/CUDAExecutor.cuh"
#include "Test.h"

using namespace Physica;
using T = float32;
using RandomSource = Random<>;
using DTensor3D = device_obj<Tensor3D<T>>;

namespace {
    void fiber() {
        const auto x = Tensor3D<T>::random_uniform<RandomSource>({4, 4, 4});
        const auto y = Tensor3D<T>::random_uniform<RandomSource>({4, 4, 4});
        const auto d_x = x.toDevice();
        const auto d_y = y.toDevice();

        const VectorND<T> answer = (x + y).fiber(1, var(), 2);
        device_obj<VectorND<T>> result = (d_x + d_y).fiber(1, var(), 2);
        expect(vectorNear(result.toHost(), answer, 1E-6));
    }

    void slice() {
        const auto x = Tensor3D<T>::random_uniform<RandomSource>({4, 4, 4});
        const auto y = Tensor3D<T>::random_uniform<RandomSource>({4, 4, 4});
        const auto d_x = x.toDevice();
        const auto d_y = y.toDevice();

        const MatrixND<T> answer = (x + y).slice(1, var(), var());
        device_obj<MatrixND<T>> result = (d_x + d_y).slice(1, var(), var());
        expect(matrixNear(result.toHost(), answer, 1E-6));
    }

    void block() {
        const auto x = Tensor3D<T>::random_uniform<RandomSource>({4, 4, 4});
        const auto y = Tensor3D<T>::random_uniform<RandomSource>({4, 4, 4});
        const auto d_x = x.toDevice();
        const auto d_y = y.toDevice();

        const Index3D from{1, 0, 2};
        const Index3D count{2, 3, 1};
        const Tensor3D<T> answer((x + y).block(from, count));
        device_obj<Tensor3D<T>> result = (d_x + d_y).block(from, count);
        expect(tensorNear(result.toHost(), answer, 1E-6));
    }

    void flatten() {
        const auto x = Tensor3D<T>::random_uniform<RandomSource>({4, 4, 4});
        const auto y = Tensor3D<T>::random_uniform<RandomSource>({4, 4, 4});
        const auto d_x = x.toDevice();
        const auto d_y = y.toDevice();

        const VectorND<T> answer = (x + y).flatten();
        device_obj<VectorND<T>> result = (d_x + d_y).flatten();
        expect(vectorNear(result.toHost(), answer, 1E-6));
    }
}

int main() {
    fiber();
    slice();
    block();
    flatten();
    return 0;
}
