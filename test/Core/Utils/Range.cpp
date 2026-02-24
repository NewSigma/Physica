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
#include <array>
#include "Physica/Core/Utils/Range.h"
#include "Physica/Core/Utils/Builtin.h"
#include "Test.h"

using namespace Physica;

namespace {
    struct array_nosize : public std::array<int, 3> {
        bool flag = true;

        [[nodiscard]] size_t size() const noexcept {
            if (flag)
                return std::array<int, 3>::size();
            unreachable();
        }
    };
}

int main() {
    std::array<int, 3> a{1, 2, 3};
    array_nosize b{};
    auto view = zip(a, b);
    using V = decltype(view);
    static_assert(std::ranges::input_range<V>);
    static_assert(std::ranges::view<V>);

    b.flag = false;
    expect(view.size() == 3);
    return 0;
}
