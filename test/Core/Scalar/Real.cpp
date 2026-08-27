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
#include "Physica/Core/Scalar/Real.h"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<>;

namespace {
    // We handle type casts in the CRTP base class. Test that we do not accidentally fall into infinite recursion.
    [[gnu::noinline]] void infiniteCompare(Scalar auto x, Scalar auto y) noexcept {
        std::ignore = x < y;
        std::ignore = x > y;
        std::ignore = x <= y;
        std::ignore = x >= y;
    }

    void real() {
        // Test x.real() returns reference
        using T = float64;
        float64 x = 0;
        x.real() += T(1);
        expect(x == T(1));
    }

    void trig() {
        using T = float64;
        constexpr uint64_t ULP = 3;
        auto x = T::random_uniform<RandomSource>();
        expect<RandomSource>(scalarNear(square(cos(x)) + square(sin(x)), T(1), ULP));
        expect<RandomSource>(scalarNear(T(1) + square(tan(x)), square(sec(x)), ULP));
        expect<RandomSource>(scalarNear(T(1) + square(cot(x)), square(csc(x)), ULP));
        expect<RandomSource>(scalarNear(square(cospi(x)) + square(sinpi(x)), T(1), ULP));
        expect<RandomSource>(scalarNear(T(1) + square(tanpi(x)), square(secpi(x)), ULP));
        expect<RandomSource>(scalarNear(T(1) + square(cotpi(x)), square(cscpi(x)), ULP));
    }
}

int main() {
    static_assert(std::formattable<float32, char>);
    static_assert(std::formattable<float64, char>);
    static_assert(std::same_as<std::strong_ordering, decltype(float32() <=> float32())>);
    infiniteCompare(float32(0), float64(0)); // Always x = y = 0, we make it a function to silent useless comparison warning
    real();
    trig();
    return 0;
}
