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
        ThreadBlock<Dynamic> block{};
        x[block.tid()] += T(block.tid());
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

    template<int NumThread>
    void syncReduce() {
        const auto d_result = device_obj<VectorND<T>>(6);
        auto kernel = [r_ = asStruct(d_result)] __device__() mutable {
            ThreadBlock<NumThread> block{};
            auto& r = r_.getDerived();
            const T index = T(block.tid());
            r[0] = block.sync_sum(index);
            r[1] = block.sync_max(index);
            r[2] = block.sync_min(index);

            const bool predicate = block.tid() < NumThread / 2;
            r[3] = T(block.sync_and(predicate));
            r[4] = T(block.sync_or(predicate));
            r[5] = T(block.sync_xor(predicate));
        };
        CUDAExecutor::launch(kernel, KernelConfig(1, NumThread));
        CUDAContext::getInstance().wait();
        const auto result = d_result.toHost();
        expect(result[0] == T(NumThread - 1) * T(NumThread) / T(2));
        expect(result[1] == T(NumThread - 1));
        expect(result[2] == T(0));
        expect(result[3] == T(0));
        expect(result[4] == (NumThread / 2 > 0 ? T(1) : T(0)));
        expect(result[5] == T(NumThread / 2 % 2));
    }
}

int main() {
    mapping1D();
    syncReduce<32>();
    // Corner case: thread count is not a power of two
    syncReduce<5>();
    syncReduce<6>();
    return 0;
}
