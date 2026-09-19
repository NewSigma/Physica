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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/DenseTensor.cuh"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Parallel/Executor/CUDAExecutor.cuh"
#include "Test.h"

using namespace Physica;
using T = float32;
using DeviceTensor = device_obj<Tensor3D<T>>;
using RandomSource = Random<>;

namespace {
    void hostDeviceCopy() {
        const auto A = Tensor3D<T>::random_uniform<RandomSource>(Index3D{4, 4, 4});
        const auto d_A = A.toDevice();
        const Tensor3D<T> B = d_A.toHost();
        expect(A.asArray() == B.asArray());
    }

    void deviceExprEval() {
        const Index3D shape{4, 4, 4};
        const Tensor3D<T> A = Tensor3D<T>::random_uniform<RandomSource>(shape);
        const Tensor3D<T> B = Tensor3D<T>::random_uniform<RandomSource>(shape);
        const Tensor3D<T> answer = A + B;

        device_obj<Tensor3D<T>> d_result = A.toDevice() + B.toDevice();
        const Tensor3D<T> result = d_result.toHost();
        expect(vectorNear(result.asArray(), answer.asArray(), 1E-6));
    }
}

int main() {
    hostDeviceCopy();
    deviceExprEval();
    return 0;
}
