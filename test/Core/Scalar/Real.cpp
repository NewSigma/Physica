/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include <limits>
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

    void constant_evaluation() {
        constexpr float64 x = 1.0;
        constexpr float64 y = 2.5;
        static_assert(x.toMachine() == 1.0);
        static_assert((x + y).toMachine() == 3.5);
        static_assert((y / x).toMachine() == 2.5);
        static_assert(static_cast<double>(y) == 2.5);
        static_assert(static_cast<float>(y) == 2.5F);
        static_assert((-y).toMachine() == -2.5);

        constexpr float32 z = 1.5F;
        static_assert(z.toMachine() == 1.5F);
        static_assert(static_cast<double>(z) == 1.5);
        static_assert(static_cast<float>(z) == 1.5F);
        static_assert((-z).toMachine() == -1.5F);
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
        constexpr double Prec = 1E-15;
        auto x = T::random_uniform<RandomSource>();
        expect<RandomSource>(scalarNear(square(cos(x)) + square(sin(x)), T(1), Prec));
        expect<RandomSource>(scalarNear(T(1) + square(tan(x)), square(sec(x)), Prec));
        expect<RandomSource>(scalarNear(T(1) + square(cot(x)), square(csc(x)), Prec));
        expect<RandomSource>(scalarNear(square(cospi(x)) + square(sinpi(x)), T(1), Prec));
        expect<RandomSource>(scalarNear(T(1) + square(tanpi(x)), square(secpi(x)), Prec));
        expect<RandomSource>(scalarNear(T(1) + square(cotpi(x)), square(cscpi(x)), Prec));
    }

    template<FloatPrec Prec>
    void subnormal() {
        using T = Real<Prec>;
        expect(T(std::numeric_limits<T>::denorm_min()).isSubNormal());
        expect(T(0).isSubNormal());
        expect(!T(1).isSubNormal());
        expect(!T(-1).isSubNormal());
        expect(!T(std::numeric_limits<T>::infinity()).isSubNormal());
        expect(!T(-std::numeric_limits<T>::infinity()).isSubNormal());
        expect(!T::nan().isSubNormal());
    }

    template<FloatPrec Prec>
    void repinf() {
        using T = Real<Prec>;
        const T inf = std::numeric_limits<T>::infinity();
        expect(reciprocal(inf).isZero());
        expect(reciprocal(-inf).isZero());
    }
}

int main() {
    static_assert(std::formattable<float32, char>);
    static_assert(std::formattable<float64, char>);
    static_assert(std::same_as<std::strong_ordering, decltype(float32() <=> float32())>);
    infiniteCompare(float32(0), float64(0)); // Always x = y = 0, we make it a function to silent useless comparison warning
    constant_evaluation();
    real();
    trig();
    subnormal<Float32>();
    subnormal<Float64>();
    repinf<Float32>();
    repinf<Float64>();
    return 0;
}
