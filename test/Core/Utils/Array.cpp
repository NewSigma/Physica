/*
 * Copyright 2025 Weibo He.
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
#include "Physica/Core/Utils/Container/Array.h"

using namespace Physica;

namespace {
template<class C>
    consteval void rangeTest() noexcept {
        using I = PtrIteratorF<C>;
        static_assert(std::indirectly_readable<I>);
        static_assert(std::indirectly_writable<I, long>);
        static_assert(std::incrementable<I>);
        static_assert(std::sized_sentinel_for<I, I>);
        static_assert(std::contiguous_iterator<I>);

        static_assert(std::ranges::sized_range<C>);
        static_assert(std::ranges::contiguous_range<C>);
        static_assert(std::ranges::common_range<C>);
        static_assert(std::ranges::viewable_range<C>);
    }
}

int main() {
    rangeTest<Array<long, 3>>();
    rangeTest<Array<long>>();
    // Test structured binding
    Array<long, 3> arr{1, 2, 3};
    auto [x, y, z] = arr;
    if (x != 1 || y != 2 || z != 3)
        return 1;
    return 0;
}
