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
        VectorND<T> result{};
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
