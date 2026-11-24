/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Scalar/Complex.h"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<>;

namespace {
    template<Scalar T, int Size>
    void range_like() {
        syntax_only([]() {
            using Pack = SIMD<T, Size>;
            static_assert(std::ranges::sized_range<Pack>);
            static_assert(std::ranges::common_range<Pack>);
            static_assert(std::ranges::random_access_range<Pack>);
            static_assert(std::random_access_iterator<std::ranges::iterator_t<Pack>>);
        });
    }

    template<Scalar T, int Size>
    void testMath(const SIMD<T, Size> x) {
        constexpr double prec = T::Prec == Float32 ? 1E-6 : 1E-15;
        /* Divide */ {
            const auto y = SIMD<T, Size>::template random_uniform<RandomSource>() + SIMD<T, Size>(1);
            const auto result = x / y;
            for (int i = 0; i < Size; ++i)
                expect<RandomSource>(scalarNear(result[i], x[i] / y[i], prec));
        }
        /* Square */ {
            const auto result = square(x);
            for (int i = 0; i < Size; ++i)
                expect<RandomSource>(scalarNear(result[i], square(x[i]), prec));
        }
        /* Exp */ {
            const auto result = exp(x);
            for (int i = 0; i < Size; ++i)
                expect<RandomSource>(scalarNear(result[i], exp(x[i]), prec));
        }
        /* ExpM1 */ {
            const auto result = expm1(x);
            for (int i = 0; i < Size; ++i)
                expect<RandomSource>(scalarNear(result[i], expm1(x[i]), prec));
        }
        /* Ln */ {
            const auto result = ln(x);
            for (int i = 0; i < Size; ++i)
                expect<RandomSource>(scalarNear(result[i], ln(x[i]), prec));
        }
        /* LnCosh */ {
            const auto result = lncosh(x);
            for (int i = 0; i < Size; ++i)
                expect<RandomSource>(scalarNear(result[i], lncosh(x[i]), prec));
        }
        if constexpr (!T::isComplex()) {
            /* Cosh */ {
                const auto result = cosh(x);
                for (int i = 0; i < Size; ++i)
                    expect(scalarNear(result[i], cosh(x[i]), prec));
            }
            /* Sinh */ {
                const auto result = sinh(x);
                for (int i = 0; i < Size; ++i)
                    expect(scalarNear(result[i], sinh(x[i]), prec));
            }
            /* Sech */
            {
                const auto result = sech(x);
                for (int i = 0; i < Size; ++i)
                    expect(scalarNear(result[i] * cosh(x[i]), T(1), prec));
            }
            /* ArcSinh */ {
                const auto result = arcsinh(x);
                for (int i = 0; i < Size; ++i)
                    expect(scalarNear(sinh(result[i]), x[i], prec));
            }
            /* ArcTanh */ {
                const auto arg = x * T(0.5);
                const auto result = arctanh(arg);
                for (int i = 0; i < Size; ++i)
                    expect(scalarNear(tanh(result[i]), arg[i], prec));
            }
            /* Sigmoid */ {
                const auto arg = x * T(2) - SIMD<T, Size>(T(0.5));
                const auto result = sigmoid(arg);
                for (int i = 0; i < Size; ++i)
                    expect(scalarNear(result[i], sigmoid(arg[i]), prec));
            }
        }
    }

    template<Scalar T, int Size>
    void test() {
        constexpr double prec = T::Prec == Float32 ? 1E-6 : 1E-15;
        auto x = SIMD<T, Size>(std::numeric_limits<T>::min());
        /* Divide */ {
            const auto result = SIMD<T, Size>(0) / x;
            for (int i = 0; i < Size; ++i)
                expect(result[i].isZero());
        }
        testMath<T, Size>(x);

        x = SIMD<T, Size>::template random_uniform<RandomSource>();
        testMath<T, Size>(x);
        /* Sqrt */ {
            const auto result = sqrt(x);
            for (int i = 0; i < Size; ++i)
                expect<RandomSource>(scalarNear(result[i], sqrt(x[i]), prec));
        }
    }
}

int main() {
    static_assert(std::formattable<SIMD<float64, 2>, char>);
    range_like<float32, 4>();
    range_like<cfloat32, 2>();

    test<float32, 4>();
    test<float32, 8>();
    test<float64, 2>();
    test<float64, 4>();
    test<cfloat32, 2>();
    test<cfloat32, 4>();
    test<cfloat64, 1>();
    test<cfloat64, 2>();
    return 0;
}
