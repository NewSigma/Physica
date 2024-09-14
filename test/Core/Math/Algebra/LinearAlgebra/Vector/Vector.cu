/*
 * Copyright 2022 Weibo He.
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
#include <Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.cuh>
#include <Physica/Core/Parallel/Executor/CUDAExecutor.cuh>

using namespace Physica;
using namespace Physica::Core;
using ScalarType = float32;
using VectorType = Vector<ScalarType>;
using DeviceVector = Core::device_obj<VectorType>;

__global__ void test_kernel(
        PlainStruct<DeviceVector> a,
        PlainStruct<DeviceVector> b,
        PlainStruct<DeviceVector> result,
        ScalarType factor) {
    result.getDerived() = a.getDerived() + b.getDerived() * factor;
}

int main() {
    const VectorType a{1, 2, 3, 4};
    auto d_a = a.toDevice();
    {
        DeviceVector d_b(4);
        d_b = d_a;
        CUDAContext::getInstance().wait();
        const VectorType b = d_b.toHost();
        if (a != b) {
            std::cout << "[Error]: Copy failed\n";
            return 1;
        }
    }
    {
        const VectorType answer = reciprocal(a);
        d_a = reciprocal(d_a);
        VectorType result;
        CUDAExecutor::wait();
        d_a.toHost(result);
        if (!vectorNear(result, answer, 1E-15)) {
            std::cout << "[Error]: Reciprocal failed\n";
            return 1;
        }
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
        DeviceVector d_result(len);
        test_kernel<<<1, len, 0, CUDAContext::getInstance()>>>(asStruct(d_a), asStruct(d_b), asStruct(d_result), factor);
        check(cudaGetLastError());
        d_result.toHostAsync(result);
        CUDAExecutor::wait();
        if (!vectorNear(result, answer, 1E-7)) {
            std::cout << "[Error]: Kernel failed\n";
            return 1;
        }
    }
    return 0;
}
