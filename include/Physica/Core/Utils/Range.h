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
#include "Physica/Core/Utils/MetaProgramming.h"

namespace Physica {
    namespace Internal {
        template<class Tuple, std::ranges::view... Vs>
        class common_zip_iterator {
            static_assert(instanceof<std::tuple, Tuple>, "[Error]: Unexpected type");
            using This = common_zip_iterator<Tuple, Vs...>;
            using StdIterator = std::ranges::iterator_t<std::ranges::zip_view<Vs...>>;
        public:
            using iterator_concept = StdIterator::iterator_concept;
            using value_type = StdIterator::value_type;
            using difference_type = StdIterator::difference_type;
        private:
            Tuple current;
        public:
            common_zip_iterator() = default;
            common_zip_iterator(Tuple current_) noexcept;
            common_zip_iterator(const This&) = default;
            common_zip_iterator(This&&) noexcept = default;
            ~common_zip_iterator() = default;
            /* Operators */
            This& operator=(const This&) = default;
            This& operator=(This&&) noexcept = default;
            constexpr This& operator++() noexcept;
            constexpr This& operator--() noexcept;
            constexpr This& operator+=(difference_type n) noexcept;
            constexpr This& operator-=(difference_type n) noexcept;
            [[nodiscard]] constexpr This operator++(int) noexcept;
            [[nodiscard]] constexpr This operator--(int) noexcept;
            [[nodiscard]] constexpr auto operator*() const noexcept;
            [[nodiscard]] constexpr auto operator[](difference_type n) const noexcept;
            [[nodiscard]] constexpr bool operator==(const This& other) const noexcept;
            [[nodiscard]] constexpr auto operator<=>(const This& other) const noexcept;
            [[nodiscard]] constexpr This operator+(difference_type n) const noexcept;
            [[nodiscard]] constexpr This operator-(difference_type n) const noexcept;
            [[nodiscard]] constexpr difference_type operator-(const This& other) const noexcept;
            /* Getters */
            template<size_t N>
            [[nodiscard]] constexpr auto& get(this auto&&) noexcept;
            /* Friends */
            friend constexpr This operator+(difference_type n, const This& ite) noexcept { return ite + n; }
            friend constexpr auto iter_move(const This& ite) noexcept {
                return std::apply([](auto&... ites) {
                    return std::make_tuple(std::ranges::iter_move(ites)...);
                }, ite.current);
            }
            friend constexpr void iter_swap(const This& lhs, const This& rhs) noexcept {
                [&]<size_t... Is>(std::index_sequence<Is...>) {
                    (std::ranges::iter_swap(std::get<Is>(lhs.current), std::get<Is>(rhs.current)), ...);
                }(std::make_index_sequence<sizeof...(Vs)>{});
            }
        };

        template<class Tuple, std::ranges::view... Vs>
        common_zip_iterator<Tuple, Vs...>::common_zip_iterator(Tuple current_) noexcept : current(std::move(current_)) {}

        template<class Tuple, std::ranges::view... Vs>
        constexpr auto common_zip_iterator<Tuple, Vs...>::operator++() noexcept -> This& {
            std::apply([](auto&... ites) {
                (++ites, ...);
            }, current);
            return *this;
        }

        template<class Tuple, std::ranges::view... Vs>
        constexpr auto common_zip_iterator<Tuple, Vs...>::operator--() noexcept -> This& {
            static_assert(std::bidirectional_iterator<StdIterator>, "[Error]: Unavailable");
            std::apply([](auto&... ites) {
                (--ites, ...);
            }, current);
            return *this;
        }

        template<class Tuple, std::ranges::view... Vs>
        constexpr auto common_zip_iterator<Tuple, Vs...>::operator+=(difference_type n) noexcept -> This& {
            static_assert(std::random_access_iterator<StdIterator>, "[Error]: Unavailable");
            std::apply([n](auto&... ites) {
                ((ites += n), ...);
            }, current);
            return *this;
        }

        template<class Tuple, std::ranges::view... Vs>
        constexpr auto common_zip_iterator<Tuple, Vs...>::operator-=(difference_type n) noexcept -> This& {
            static_assert(std::random_access_iterator<StdIterator>, "[Error]: Unavailable");
            std::apply([n](auto&... ites) {
                ((ites -= n), ...);
            }, current);
            return *this;
        }

        template<class Tuple, std::ranges::view... Vs>
        constexpr auto common_zip_iterator<Tuple, Vs...>::operator++(int) noexcept -> This {
            static_assert(std::forward_iterator<StdIterator>, "[Error]: Unavailable");
            auto tmp = *this;
            ++(*this);
            return tmp;
        }

        template<class Tuple, std::ranges::view... Vs>
        constexpr auto common_zip_iterator<Tuple, Vs...>::operator--(int) noexcept -> This {
            static_assert(std::bidirectional_iterator<StdIterator>, "[Error]: Unavailable");
            auto tmp = *this;
            --(*this);
            return tmp;
        }

        template<class Tuple, std::ranges::view... Vs>
        constexpr auto common_zip_iterator<Tuple, Vs...>::operator*() const noexcept {
            return std::apply([](auto&... ites) {
                return std::make_tuple(*ites...);
            }, current);
        }

        template<class Tuple, std::ranges::view... Vs>
        constexpr auto common_zip_iterator<Tuple, Vs...>::operator[](difference_type n) const noexcept {
            static_assert(std::random_access_iterator<StdIterator>, "[Error]: Unavailable");
            return std::apply([n](auto&... ites) -> decltype(auto) {
                return std::make_tuple(ites[n]...);
            }, current);
        }

        template<class Tuple, std::ranges::view... Vs>
        constexpr bool common_zip_iterator<Tuple, Vs...>::operator==(const This& other) const noexcept {
            return std::get<0>(current) == std::get<0>(other.current);
        }

        template<class Tuple, std::ranges::view... Vs>
        constexpr auto common_zip_iterator<Tuple, Vs...>::operator<=>(const This& other) const noexcept {
            return std::get<0>(current) <=> std::get<0>(other.current);
        }

        template<class Tuple, std::ranges::view... Vs>
        constexpr auto common_zip_iterator<Tuple, Vs...>::operator+(difference_type n) const noexcept -> This {
            static_assert(std::random_access_iterator<StdIterator>, "[Error]: Unavailable");
            auto r = *this;
            r += n;
            return r;
        }

        template<class Tuple, std::ranges::view... Vs>
        constexpr auto common_zip_iterator<Tuple, Vs...>::operator-(difference_type n) const noexcept -> This {
            static_assert(std::random_access_iterator<StdIterator>, "[Error]: Unavailable");
            auto r = *this;
            r -= n;
            return r;
        }

        template<class Tuple, std::ranges::view... Vs>
        constexpr auto common_zip_iterator<Tuple, Vs...>::operator-(const This& other) const noexcept -> difference_type {
            return std::get<0>(current) - std::get<0>(other.current);
        }

        template<class Tuple, std::ranges::view... Vs>
        template<size_t N>
        constexpr auto& common_zip_iterator<Tuple, Vs...>::get(this auto&& self) noexcept {
            return std::get<N>(self.current);
        }

        template<std::ranges::view... Vs>
        class common_zip final : public std::ranges::view_interface<common_zip<Vs...>> {
            static_assert((std::ranges::input_range<Vs> && ...));
            static_assert((std::ranges::common_range<Vs> && ...), "[Error]: common_zip requires common_range");
            static_assert(sizeof...(Vs) > 1, "[Error]: Unnecessary call to zip");
            using This = common_zip<Vs...>;
        private:
            std::tuple<Vs...> views;
        public:
            constexpr explicit common_zip(Vs... views_) noexcept;
            common_zip(const This&) = default;
            common_zip(This&&) noexcept = default;
            ~common_zip() = default;
            /* Operators */
            This& operator=(const This&) = default;
            This& operator=(This&&) noexcept = default;
            /* Operations */
            [[nodiscard]] constexpr auto begin(this auto&&) noexcept;
            [[nodiscard]] constexpr auto end(this auto&&) noexcept;
            [[nodiscard]] constexpr auto size() const noexcept;
            /* Getters */
            template<size_t N>
            [[nodiscard]] constexpr auto& get(this auto&&) noexcept;
        };

        template<std::ranges::view... Vs>
        constexpr auto common_zip<Vs...>::begin(this auto&& self) noexcept {
            auto current = std::apply([](auto&... views) {
                return std::make_tuple(std::ranges::begin(std::forward<decltype(views)>(views))...);
            }, self.views);
            return common_zip_iterator<decltype(current), Vs...>(std::move(current));
        }

        template<std::ranges::view... Vs>
        constexpr auto common_zip<Vs...>::end(this auto&& self) noexcept {
            if constexpr ((std::ranges::random_access_range<Vs> && ...))
                return self.begin() + self.size();
            else {
                auto current = std::apply([](auto&... views) {
                    return std::make_tuple(std::ranges::end(views)...);
                }, self.views);
                return common_zip_iterator<decltype(current), Vs...>(std::move(current));
            }
        }

        template<std::ranges::view... Vs>
        constexpr common_zip<Vs...>::common_zip(Vs... views_) noexcept : views(std::move(views_)...) {
            std::apply([](const auto& view0, const auto&... rest) {
                auto size = std::ranges::size(view0);
                assert(((std::ranges::size(rest) == size) && ...) && "[Error]: Sizes mismatch");
                static_assert(((std::same_as<decltype(size), decltype(std::ranges::size(rest))>) && ...), "[Error]: Size type mismatch");
            }, views);
        }

        template<std::ranges::view... Vs>
        constexpr auto common_zip<Vs...>::size() const noexcept {
            static_assert((std::ranges::sized_range<Vs> && ...), "[Error]: Unavailable");
            return std::get<0>(views).size();
        }

        template<std::ranges::view... Vs>
        template<size_t N>
        constexpr auto& common_zip<Vs...>::get(this auto&& self) noexcept {
            return std::get<N>(self.views);
        }
    }
    /**
     * 1. Assume the views have identical size here and avoid the std::min call in std::ranges::views::zip_view.
     * 2. Allow unzip into indivisual views/iterators
     */
    [[nodiscard]] constexpr auto zip(std::ranges::viewable_range auto&&... ranges) noexcept {
        return Internal::common_zip<std::ranges::views::all_t<decltype(ranges)>...>(std::forward<decltype(ranges)>(ranges)...);
    }
}

namespace std {
    template<ranges::view... Vs>
    struct tuple_size<Physica::Internal::common_zip<Vs...>> : public integral_constant<size_t, sizeof...(Vs)> {};

    template<size_t N, ranges::view... Vs>
    struct tuple_element<N, Physica::Internal::common_zip<Vs...>> {
        using type = tuple_element_t<N, tuple<Vs...>>;
    };

    template<class Tuple, ranges::view... Vs>
    struct tuple_size<Physica::Internal::common_zip_iterator<Tuple, Vs...>> : public integral_constant<size_t, sizeof...(Vs)> {};

    template<size_t N, class Tuple, ranges::view... Vs>
    struct tuple_element<N, Physica::Internal::common_zip_iterator<Tuple, Vs...>> {
        using type = tuple_element_t<N, Tuple>;
    };
}
