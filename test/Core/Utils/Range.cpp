/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <array>
#include <concepts>
#include <utility>
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
    array_nosize b{{3, 1, 2}};
    auto view = zip(a, b);
    /* Test size */ {
        using V = decltype(view);
        static_assert(std::same_as<std::ranges::range_value_t<V>, std::pair<int, int>>, "[Error]: Prefer std::pair to tuple2D, std::pair generally provides better codegen");
        static_assert(std::ranges::input_range<V>);
        static_assert(std::ranges::view<V>);

        b.flag = false;
        expect(view.size() == 3);
    }
    /* unzip */ {
        auto [va, vb] = view;
        expect(va[2] == vb[0]);

        auto [ia, ib] = ++view.begin();
        expect(*ia == 2 && *ib == 1);
    }
    return 0;
}
