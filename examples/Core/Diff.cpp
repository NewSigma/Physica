/*
 * Copyright 2025-2026 Weibo He.
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
#include <print>
#include "Physica/Core/Scalar/Diff.h"

using namespace Physica;

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
}
/**
 * We demonstrate basic concepts of automatic differentiation implemented in Physica.
 */
int main() {
    std::println("{} = {}", forward(), reverse());
    return 0;
}
