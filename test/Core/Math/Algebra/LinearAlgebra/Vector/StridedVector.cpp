/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DiffVector.h"
#include "Physica/Core/Scalar/Complex.h"
#include "Test.h"

using namespace Physica;

int main() {
    syntax_only([]() {
        using R = decltype(VectorND<cfloat64>{}.reals().view());
        static_assert(std::ranges::sized_range<R>);
        static_assert(std::ranges::common_range<R>);
        static_assert(std::ranges::random_access_range<R>);
        static_assert(std::ranges::view<R>);
    });
    syntax_only([]() {
        using dcfloat64 = Diff<cfloat64, DiffMode::Forward, 2>;
        VectorND<dcfloat64> v(2);
        std::ignore = v.reals()[0];
        std::ignore = v.imags()[0];

        VectorND<cfloat64> v1(2);
        std::ignore = v1.imags()[0];
    });
    return 0;
}
