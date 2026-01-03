/*
 * Copyright 2023-2025 Weibo He.
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
#include "Physica/Core/Scalar/Diff.h"
#include "Test.h"

using namespace Physica;
using T = float64;

namespace {
    void testFunc() {
        bool good = true;
        {
            using dfloat = Diff<T, DiffMode::Forward, 1>;
            auto func = [](dfloat x, dfloat y) -> dfloat {
                return square(x - T(1.0)) + square(y - T(2.0));
            };
            const T x = 3;
            const T y = 4;
            const dfloat result = func(dfloat(x, 1), dfloat(y, 1));
            const T answer = (x + y - 3.0) * 2.0;
            good &= scalarNear(result.grad(), answer, 1E-15);
        }
        {
            using dfloat = Diff<T, DiffMode::Forward, 2>;
            dfloat x{3, 1};
            dfloat y = square(x);
            good &= scalarNear(y.template grad<2>(), float64(2), 1E-15);
        }
        expect(good);
    }

    void testMath() {
        bool good = true;
        {
            using dfloat = Diff<T, DiffMode::Forward, 1>;
            auto result = T(3) / dfloat(T(2), T(4));
            good &= scalarNear(result.value(), T(1.5), 1E-15);
            good &= scalarNear(result.grad(), T(-3), 1E-15);
        }
        {
            using dfloat = Diff<T, DiffMode::Forward, 2>;
            dfloat x(3, 1);
            auto y = reciprocal(x);
            good &= scalarNear(y.grad().value(), -square(reciprocal(x.value())), 1E-15);
            good &= scalarNear(y.grad<2>(), pow(reciprocal(x.value()), T(3)) * T(2), 1E-15);

            y = sqrt(x);
            good &= scalarNear(y.grad().value(), reciprocal(T(2) * sqrt(x.value())), 1E-15);
            good &= scalarNear(y.grad<2>(), -reciprocal(T(4) * x.value() * sqrt(x.value())), 1E-15);
        }
        expect(good);
    }

    void testSIMD() {
        T value[4]{1.5, -1.5, 0, 2};
        T grad1[4]{1, 1, 1, 0};
        T grad2[4]{0, 0, 0, 0};

        using dfloat = Diff<T, DiffMode::Forward, 2>;
        SIMD<dfloat, 4> packet{};
        packet.load({value, {grad1, grad2}});

        bool good = true;
        auto result = abs(packet);
        for (int i = 0; i < 4; ++i)
            good &= scalarNear(abs(dfloat(value[i], grad1[i])), result[i], 1E-15);

        result = square(packet);
        for (int i = 0; i < 4; ++i)
            good &= scalarNear(square(dfloat(value[i], grad1[i])), result[i], 1E-15);

        result = reciprocal(packet);
        for (int i = 0; i < 4; ++i) {
            if (value[i].isZero())
                continue;
            good &= scalarNear(reciprocal(dfloat(value[i], grad1[i])), result[i], 1E-15);
        }

        result = exp(packet);
        for (int i = 0; i < 4; ++i)
            good &= scalarNear(exp(dfloat(value[i], grad1[i])), result[i], 1E-15);

        expect(good);
    }

    void testReverse() {
        using dfloat = Diff<float64, DiffMode::Reverse>;
        auto x = dfloat(2);
        /* Simple */ {
            const auto y = sin(x).reverse();
            expect(y == sin(x.value()));
            expect(x.grad() == cos(x.value()));
        }
        /* Test r-value 1 */ {
            x.zero_grad();
            const auto y = sin(sin(x)).reverse();
            expect(y.value() == sin(sin(x.value())));
            expect(x.grad() == cos(x.value()) * cos(sin(x.value())));
        }
        /* Test r-value 2 */ {
            auto func = [](const dfloat& x) {
                return sin(sin(x));
            };
            x.zero_grad();
            func(x).reverse();
            expect(x.grad() == cos(x.value()) * cos(sin(x.value())));
        }
        {
            auto func = [](dfloat& x, dfloat& y) {
                return square(x - T(1.0)) + square(y - T(2.0));
            };
            dfloat x(3);
            dfloat y(4);
            func(x, y).reverse();
            expect(scalarNear(x.grad(), (x.value() - 1.0) * 2.0, 1E-15));
            expect(scalarNear(y.grad(), (y.value() - 2.0) * 2.0, 1E-15));
        }
    }
}

int main() {
    static_assert(std::formattable<Diff<float64, DiffMode::Forward>, char>);
    static_assert(std::formattable<DiffCoro<Diff<float64, DiffMode::Reverse>>, char>);
    testFunc();
    testMath();
    testSIMD();
    testReverse();
    return 0;
}
