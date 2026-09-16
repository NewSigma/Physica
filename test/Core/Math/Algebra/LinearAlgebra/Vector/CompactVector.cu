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
#include "Physica/Core/Parallel/Executor/CUDAExecutor.cuh"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<>;

namespace {
    void read() {
        using Tr = float32;
        using Tc = cfloat32;
        const auto real = VectorND<Tr>::random_uniform<RandomSource>(12);
        auto& ctx = CUDAContext::getInstance();

        device_obj<VectorND<Tc>> d_vec(6);
        d_vec.read(real);
        ctx.wait();

        VectorND<Tr> back(real.getLength(), 0);
        back.read(d_vec);
        ctx.wait();
        expect(back == real);
    }
}

int main() {
    read();
    return 0;
}
