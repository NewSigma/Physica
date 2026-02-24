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
#pragma once

#include <cassert>
#include <ranges>

namespace Physica {
    /**
     * We assume the views have identical size here and avoid the std::min call in std::ranges::views::zip_view.
     */
    constexpr auto zip(std::ranges::viewable_range auto&&... ranges) noexcept {
        static_assert(sizeof...(ranges) > 1, "[Error]: Unnecessary call to zip");
        [](const auto& view0, const auto&... rest) {
            auto size = view0.size();
            assert(((rest.size() == size) && ...) && "[Error]: Sizes mismatch");
            static_assert(((std::same_as<decltype(size), decltype(rest.size())>) && ...), "[Error]: Size type mismatch");
        }(ranges...);

        using Base = std::ranges::zip_view<std::ranges::views::all_t<decltype(ranges)>...>;
        class common_zip : public Base {
        public:
            using Base::Base;

            constexpr auto size() const noexcept {
                return std::get<0>(reinterpret_cast<const std::tuple<std::ranges::views::all_t<decltype(ranges)>...>&>(*this)).size();
            }
        };
        return common_zip(std::forward<decltype(ranges)>(ranges)...);
    }
}
