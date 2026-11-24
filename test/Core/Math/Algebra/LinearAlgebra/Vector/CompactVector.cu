/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
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
