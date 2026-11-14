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
#include "Physica/Core/Math/Algebra/Canonicalization.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Test.h"

using namespace Physica;

using Expr1 = decltype(VectorND<float64>() * float64(1));
using Expr2 = decltype(VectorND<float64>() + VectorND<float64>());
using Expr3 = decltype(VectorND<float32>() + VectorND<float32>());

namespace {
    template<class T1, class T2>
    consteval bool exclusive() noexcept {
        return canonicalized<T1, T2>() != canonicalized<T2, T1>();
    }

    template<class T1, class T2, class T3>
    consteval bool transitive() noexcept {
        return !(canonicalized<T1, T2>() && canonicalized<T2, T3>()) || canonicalized<T1, T3>();
    }
}

static_assert(Internal::type_order<VectorND<float64>, VectorND<float64>>() == 0);
static_assert(Internal::type_order<VectorND<float64>, VectorND<float32>>() != 0);
static_assert((Internal::type_order<int, long>() < 0) && (Internal::type_order<long, int>() > 0));
// Scalars with equal scores
static_assert(exclusive<float64, float32>());
// Vectors with equal scores
static_assert(exclusive<VectorND<float64>, VectorND<float32>>());
// Equal scores are broken by the type order
static_assert(canonicalized<VectorND<float64>, VectorND<float32>>() == (Internal::type_order<VectorND<float64>, VectorND<float32>>() > 0));
// Expression templates with equal scores
static_assert(exclusive<Expr1, Expr2>());
static_assert(exclusive<Expr2, Expr3>());
// Mixed concept scores are ordered by concept score
static_assert(canonicalized<VectorND<float64>, float64>());
static_assert(canonicalized<Expr1, VectorND<float64>>());
static_assert(transitive<float64, VectorND<float64>, VectorND<float32>>());
static_assert(transitive<VectorND<float32>, VectorND<float64>, Expr1>());
static_assert(transitive<VectorND<float64>, Expr2, Expr3>());

int main() {
    // Reordering operands must yield the same expression type
    const auto v64 = VectorND<float64>();
    const auto v32 = VectorND<float32>();
    static_assert(std::is_same_v<decltype(hadamard(v64, v32)), decltype(hadamard(v32, v64))>);
    static_assert(std::is_same_v<decltype(hadamard(v64, v64)), decltype(hadamard(v64, v64))>);

    syntax_only([]() {
        auto expr = VectorND<float64>() * float64(1);
        auto complex_vec = VectorND<cfloat64>();
        std::ignore = hadamard(expr, complex_vec);
    });
}
