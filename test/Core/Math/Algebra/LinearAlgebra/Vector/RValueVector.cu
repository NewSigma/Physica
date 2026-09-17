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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DiffVector.cuh"
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

    void transpose_hermite() {
        expect(device_obj<Vector1D<float32>>{}.transpose().getOrder() == 1);
        expect(device_obj<Vector1D<cfloat32>>{}.hermite().getOrder() == 1);
    }

    void equality() {
        {
            VectorND<float64> a = {1, 2, 3, 4};
            VectorND<float64> expected = {2, 4, 6, 8};
            const auto d_a = a.toDevice();
            const auto d_expected = expected.toDevice();
            const auto expr = d_a * float64(2);
            expect(expr == d_expected);
            expect(d_a != expr);
        }
        {
            const auto big = VectorND<float64>::random_uniform<RandomSource>(1024 + 1);
            const auto d_big = big.toDevice();
            const auto copy = d_big;
            expect(d_big == copy);
            expect(d_big != copy * float64(2));
        }
        using dfloat = Diff<float64, DiffMode::Forward, 1>;
        using DiffDevice = device_obj<VectorND<dfloat>>;
        const DiffDevice d_a(3, float64(0));
        DiffDevice d_b(3, float64(0));
        expect(d_a == d_b);
        expect(d_a.values() == d_b.values());

        d_b.grads() = float64(1);
        expect(d_a != d_b);
        expect(d_a.values() == d_b.values());
    }
}

int main() {
    host_dev_dot();
    cooperative_dot();
    transpose_hermite();
    equality();
    return 0;
}
