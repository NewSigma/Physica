/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<>;

namespace {
    void assign_range() noexcept {
        using T = float32;
        auto pred = [](float32 x) static noexcept { return x != T(5); };
        auto x = VectorND<T>{1, 2, 3, 4, 5, 6, 7, 8};
        auto v = x.view() | std::views::filter(pred);
        VectorND<T> result = VectorND<T>::junk(x.getLength());
        result = v;
        expect(result.getLength() == x.getLength() - 1);
        for (T elem : result)
            expect(pred(elem));
    }
}

int main() {
    assign_range();
    return 0;
}
