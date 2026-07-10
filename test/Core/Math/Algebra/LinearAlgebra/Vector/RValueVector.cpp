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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/IdentityMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DiffVector.h"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<>;

namespace {
    void sum() {
        using T = float32;
        const auto x = VectorND<T>::random_uniform<Random<>>(16);
        expect<RandomSource>(scalarNear(x.sum(), std::ranges::fold_left(x.view(), T(0), std::plus<>()), 3UL));

        const auto y = VectorND<T>::random_uniform<Random<>>(17);
        expect<RandomSource>(scalarNear(y.sum(), std::ranges::fold_left(y.view(), T(0), std::plus<>()), 3UL));
    }

    void argminTest() noexcept {
        const auto x = VectorND<float32>::random_uniform<Random<>>(16);
        expect<RandomSource>(x.min() == x[x.argmin()]);
    }

    void reshapeTest() {
        using T = float32;
        auto id = IdentityMatrix<T>(3);
        auto x = VectorND<T>::random_uniform<RandomSource>(id.getSize());
        MatrixND<T> y = id * (x + x).reshape_col(3, 3);
        expect<RandomSource>(vectorNear(y.flatten(), x * T(2), uint64_t(0)));

        // Test that we fallback to Col major if cannot infer major from input matrix
        syntax_only([]() {
            auto mat = IdentityMatrix<T>{};
            auto re = abs(VectorND<T>{}).reshape_like(mat);
            static_assert(mat.isBothMajor() && re.isColMatrix());
        });
    }

    void dot_reverse() {
        using T = float64;
        using dfloat = Diff<T, DiffMode::Reverse>;
        VectorND<dfloat> x{1, 2, 3};
        {
            auto result = (square(x) * VectorND<T>{1, 1, 1});
            result.reverse();
        }
        expect(x.grads() == T(2) * x.values());
    }

    void transpose_hermite() {
        static_assert(DenseVector<float64>{}.transpose().getRowAtCompile() == 1);
        static_assert(DenseVector<cfloat64>{}.hermite().getRowAtCompile() == 1);
        expect(abs(Vector1D<float64>{}).transpose().getOrder() == 1);
        expect(abs(Vector1D<cfloat64>{}).hermite().getOrder() == 1);
    }
}

int main() {
    syntax_only([]() {
        VectorND<float32> a, b;
        using R = decltype((a + b).view());
        static_assert(std::ranges::sized_range<R>);
        static_assert(std::ranges::common_range<R>);
        static_assert(std::ranges::random_access_range<R>);
        static_assert(std::ranges::constant_range<R>);
        static_assert(std::ranges::view<R>);
    });

    sum();
    argminTest();
    reshapeTest();
    dot_reverse();
    transpose_hermite();
    return 0;
}
