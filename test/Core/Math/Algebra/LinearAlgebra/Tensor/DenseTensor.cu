/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/DenseTensor.cuh"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Parallel/Executor/CUDAExecutor.cuh"
#include "Test.h"

using namespace Physica;
using T = float32;
using DTensor3D = device_obj<Tensor3D<T>>;
using RandomSource = Random<>;

namespace {
    void hostDeviceCopy() {
        const auto A = Tensor3D<T>::random_uniform<RandomSource>(Index3D{4, 4, 4});
        const auto d_A = A.toDevice();
        const Tensor3D<T> B = d_A.toHost();
        expect(A.asArray() == B.asArray());
    }

    void deviceExprEval() {
        const Index3D shape{4, 4, 4};
        const Tensor3D<T> A = Tensor3D<T>::random_uniform<RandomSource>(shape);
        const Tensor3D<T> B = Tensor3D<T>::random_uniform<RandomSource>(shape);
        const Tensor3D<T> answer = A + B;

        device_obj<Tensor3D<T>> d_result = A.toDevice() + B.toDevice();
        const Tensor3D<T> result = d_result.toHost();
        expect(vectorNear(result.asArray(), answer.asArray(), 1E-6));
    }
}

static_assert(DTensor3D::isCompact(), "DenseTensor is a compact object");
static_assert(DTensor3D::isStrided(), "DenseTensor is a strided object");
static_assert(DTensor3D::getStrideAtCompile()[0] == Dynamic, "Dynamic shape strides are unknown at compile time");
static_assert(DTensor3D::getStrideAtCompile()[2] == 1, "DenseTensor is row-major");
static_assert(DTensor3D::getStrideAtCompile(0) == DTensor3D::getStrideAtCompile()[0]);
static_assert(DTensor3D::getStrideAtCompile(2) == DTensor3D::getStrideAtCompile()[2]);

int main() {
    hostDeviceCopy();
    deviceExprEval();
    return 0;
}
