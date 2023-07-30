/*
 * Copyright 2022 WeiBo He.
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
#include <random>
#include <iostream>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.cuh"
#include "Physica/Core/Parallel/Executor/CudaExecutor.cuh"

using namespace Physica;
using namespace Physica::Core;
using ScalarType = Scalar<Float>;
using VectorType = Vector<ScalarType>;

__global__ void test_kernel(
        PlainStruct<device_obj<VectorType>> a,
        PlainStruct<device_obj<VectorType>> b,
        PlainStruct<device_obj<VectorType>> result,
        ScalarType factor) {
    result.getDerived() = a.getDerived() + b.getDerived() * factor;
}

int main() {
    VectorType a{1, 2, 3, 4};
    auto d_a = a.toDevice();
    device_obj<VectorType> d_b(4);
    d_b = d_a;
    cudaCheck(cudaStreamSynchronize(nullptr));
    VectorType b = d_b.toHost();
    if (a != b)
        return 1;
    {
        const VectorType answer = reciprocal(a);
        d_b = reciprocal(d_b);
        cudaCheck(cudaStreamSynchronize(nullptr));
        VectorType result;
        d_b.toHost(result);
        if (!vectorNear(result, answer, 1E-15))
            return 1;
    }
    {
        constexpr size_t len = 32;
        const ScalarType factor = 1.4;
        std::mt19937 gen{};
        const auto a = VectorType::random_uniform(len, gen);
        const auto b = VectorType::random_uniform(len, gen);
        const VectorType answer = a + b * factor;
        auto d_a = a.toDevice();
        auto d_b = b.toDevice();
        VectorType result(len);
        device_obj<VectorType> d_result(len);
        test_kernel<<<1, len>>>(asStruct(d_a), asStruct(d_b), asStruct(d_result), factor);
        d_result.toHostAsync(result);
        CudaExecutor::wait();
        if (!vectorNear(result, answer, 1E-7))
            return 1;
    }
    return 0;
}
