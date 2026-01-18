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
#include "Test.h"

using namespace Physica;
using RandomSource = Random<MCG>;
using T = float32;

namespace {
    __device__ void add1D(device_obj<VectorND<T>>& x) {
        ThreadBlock block{};
        x[block.rank()] += T(block.rank());
    }
    // Test that ThreadBlock changes logical threads layout
    void mapping1D() {
        auto data = VectorND<T>::random_uniform<RandomSource>(64);
        auto d_data = data.toDeviceAsync();
        auto kernel = [x = asStruct(d_data)] __device__() mutable {
            add1D(x.getDerived());
        };
        CUDAExecutor::launch(kernel, KernelConfig(1, {32, 2}));
        CUDAContext::getInstance().wait();
        VectorND<T> data1 = d_data.toHost();
        for (size_t i = 0; i < data.getLength(); ++i)
            expect(data1[i] == data[i] + T(i));
    }
}

int main() {
    mapping1D();
    return 0;
}
