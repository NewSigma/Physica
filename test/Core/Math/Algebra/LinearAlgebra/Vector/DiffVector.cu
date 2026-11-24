/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DiffVector.cuh"
#include "Physica/Core/Math/Random/Random.h"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<>;
using T = float32;

int main() {
    using VectorType = device_obj<VectorND<Diff<T, DiffMode::Reverse>>>;
    auto v = VectorType(8, 0);
    v.reverse(T(1));
    CUDAContext::getInstance().wait();
    const auto grads = v.grads().toHost();

    for (auto elem : grads)
        expect(elem == T(1));
    return 0;
}
