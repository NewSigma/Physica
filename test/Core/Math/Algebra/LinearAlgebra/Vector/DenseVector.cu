/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.cuh"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Parallel/Executor/CUDAExecutor.cuh"
#include "Test.h"

using namespace Physica;
using ScalarType = float32;
using VectorType = VectorND<ScalarType>;
using DeviceVector = device_obj<VectorType>;
using RandomSource = Random<MT19937, std::mt19937::default_seed>;

namespace {
    void transpose_hermite() {
        expect(device_obj<Vector1D<float32>>{}.transpose().getOrder() == 1);
        expect(device_obj<Vector1D<cfloat32>>{}.hermite().getOrder() == 1);
    }
}

static_assert(device_obj<VectorND<float32>>{}.transpose().isCompact(), "Transpose of a compact vector is a compact vector");

int main() {
    const VectorType a{1, 2, 3, 4};
    auto d_a = a.toDevice();
    {
        DeviceVector d_b(4);
        d_b = d_a;
        CUDAContext::getInstance().wait();
        const VectorType b = d_b.toHost();
        expect(a == b);
    }
    {
        const VectorType answer = reciprocal(a);
        d_a = reciprocal(d_a);
        VectorType result;
        CUDAContext::getInstance().wait();
        d_a.toHost(result);
        expect(vectorNear(result, answer, 1E-15));
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
            int i = static_cast<int>(threadIdx.x);
            result.getDerived()[i] = a.getDerived()[i] + b.getDerived()[i] * factor;
        }, KernelConfig(1, len));
        CUDAExecutor::wait();

        d_result.toHost(result);
        expect(vectorNear(result, answer, 1E-6));
    }
    transpose_hermite();
    return 0;
}
