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
#include "Physica/Core/Math/Random/Random.h"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<>;

namespace {
    void host_dev_dot() {
        constexpr int NumThread = 32;
        using T = float32;
        const auto v = VectorND<T>::random_uniform<RandomSource>(NumThread);
        const auto dv = v.toDevice();
        expect<RandomSource>(scalarNear(v * v, dv * dv, 2UL));
    }

    void cooperative_dot() {
        constexpr int NumThread = 32;
        using T = float32;
        const auto v = VectorND<T>::random_uniform<RandomSource>(NumThread);
        auto d_v = v.toDevice();
        CUDAExecutor::launch([v_ = asStruct(d_v)] __device__() mutable {
            auto& v = v_.getDerived();
            T v2 = dot(v, v).calc(ThreadBlock<NumThread>{});
            v[threadIdx.x] = v2;
        }, {1, NumThread});
        expect<RandomSource>(scalarNear(v * v, d_v.toHost()[0], 2UL));
    }
}

int main() {
    host_dev_dot();
    cooperative_dot();
    return 0;
}
