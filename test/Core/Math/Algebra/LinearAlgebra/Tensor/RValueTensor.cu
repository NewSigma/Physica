/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
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

        const VectorND<T> answer = (x + y).fiber(1UZ, var(), 2UZ);
        device_obj<VectorND<T>> result = (d_x + d_y).fiber(1UZ, var(), 2UZ);
        expect(vectorNear(result.toHost(), answer, 1E-6));
    }

    void slice() {
        const auto x = Tensor3D<T>::random_uniform<RandomSource>({4, 4, 4});
        const auto y = Tensor3D<T>::random_uniform<RandomSource>({4, 4, 4});
        const auto d_x = x.toDevice();
        const auto d_y = y.toDevice();

        const MatrixND<T> answer = (x + y).slice(1UZ, var(), var());
        device_obj<MatrixND<T>> result = (d_x + d_y).slice(1UZ, var(), var());
        expect(matrixNear(result.toHost(), answer, 1E-6));
    }

    void block() {
        {
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
        {
            const auto x = Tensor4D<T>::random_uniform<RandomSource>({2, 3, 4, 5});
            const auto y = Tensor4D<T>::random_uniform<RandomSource>({2, 3, 4, 5});
            const auto d_x = x.toDevice();
            const auto d_y = y.toDevice();

            const Tensor4D<T> answer((x + y).block(1UZ, var(), var(), var()));
            device_obj<Tensor4D<T>> result = (d_x + d_y).block(1UZ, var(), var(), var());
            expect(tensorNear(result.toHost(), answer, 1E-6));
        }
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
