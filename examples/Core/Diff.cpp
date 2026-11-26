/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include <print>
#include "Physica/Core/Scalar/Diff.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiffDenseMatrix.h"

using namespace Physica;
using T = float64;
using dfloat = Diff<T, DiffMode::Reverse>;

namespace {
    float64 forward() noexcept {
        using dfloat = Diff<float64, DiffMode::Forward>;
        dfloat x(0.5, 1); // x = 0.5, dx/dx = 1
        dfloat y = sin(x) / x; // Compute value as well as grad
        return y.grad();
    }

    float64 reverse() noexcept {
        using dfloat = Diff<float64, DiffMode::Reverse>;
        dfloat x(0.5); // x = 0.5, initial grad is 0
        {
            auto y = sin(x) / x; // Compute value, push compute graph node
            static_assert(std::same_as<decltype(y), CoDiff<dfloat>>, "Struct CoDiff<dfloat> contains a dfloat and a handle to compute graph node");
            y.reverse(1); // dy/dy = 1, prepare for lazy reverse propagating
        } // Pop compute graph node, propagate grads
        return x.grad(); // x receives grad
    }

    void linalg() {
        MatrixND<dfloat> A{{1, 2}, {-1, 1}};
        VectorND<dfloat> x{1, 1};
        VectorND<dfloat> b{-1, -0.5};

        auto y = relu(A * x + b);
        auto loss = square(y).sum();
        loss.reverse_final(); // Immediate propagation without waiting for y to be destroyed

        std::println("dL/dA =\n{}", A.grads());
        std::println("dL/dx = {}", x.grads());
        std::println("dL/db = {}", b.grads());
    }
}
/**
 * We demonstrate basic concepts of automatic differentiation implemented in Physica.
 */
int main() {
    std::println("{} = {}\n", forward(), reverse());
    linalg();
    return 0;
}
