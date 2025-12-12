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
#include <iostream>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.cuh"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Parallel/Executor/CUDAExecutor.cuh"

using namespace Physica;
using ScalarType = float32;
using VectorType = VectorND<ScalarType>;
using DeviceVector = device_obj<VectorType>;
using RandomSource = Random<MT19937, std::mt19937::default_seed>;

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
        CUDAContext::getInstance().wait();
        d_a.toHost(result);
        if (!vectorNear(result, answer, 1E-15)) {
            std::cout << "[Error]: Reciprocal failed\n";
            return 1;
        }
    }
    {
        constexpr size_t len = 32;
        const ScalarType factor = 1.4;
        const auto a = VectorType::random_uniform<RandomSource>(len);
        const auto b = VectorType::random_uniform<RandomSource>(len);
        const VectorType answer = a + b * factor;
        auto d_a = a.toDevice();
        auto d_b = b.toDevice();
        VectorType result(len);
        DeviceVector d_result(len);
        CUDAExecutor::launch([a = asStruct(d_a), b = asStruct(d_b), result = asStruct(d_result), factor] __device__() mutable {
            int i = threadIdx.x;
            result.getDerived()[i] = a.getDerived()[i] + b.getDerived()[i] * factor;
        }, KernelConfig(1, len)).wait();

        d_result.toHost(result);
        if (!vectorNear(result, answer, 1E-6)) {
            std::cout << "[Error]: Bad result\n";
            return 1;
        }
    }
    return 0;
}
