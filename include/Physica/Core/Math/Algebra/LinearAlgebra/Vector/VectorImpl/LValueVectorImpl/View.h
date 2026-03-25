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

#include "../LValueVector.h"

namespace Physica {
    template<class Derived>
    template<Vector V>
    class LValueVector<Derived>::View : public std::ranges::view_base {
        static_assert(!std::is_reference_v<V>);
        using This = View;

        class Iterator;

        V* vec;
    public:
        [[gnu::always_inline]] constexpr View(V& vec) noexcept : vec(&vec) {}
        constexpr View(const This&) = default;
        constexpr View(This&&) noexcept = default;
        constexpr ~View() = default;
        /* Operators */
        constexpr This& operator=(const This&) = default;
        constexpr This& operator=(This&&) noexcept = default;
        /* Operations */
        [[nodiscard, gnu::always_inline]] constexpr auto begin(this auto&&) noexcept;
        [[nodiscard, gnu::always_inline]] constexpr auto end(this auto&&) noexcept;
        /* Getters */
        [[nodiscard, gnu::always_inline]] constexpr V& vector() const noexcept { return *vec; }
        [[nodiscard, gnu::always_inline]] constexpr size_t size() const noexcept { return vec->getLength(); }
    };

    template<class Derived>
    template<Vector V>
    constexpr auto LValueVector<Derived>::View<V>::begin(this auto&& self) noexcept {
        return Iterator(&self, 0);
    }

    template<class Derived>
    template<Vector V>
    constexpr auto LValueVector<Derived>::View<V>::end(this auto&& self) noexcept {
        return Iterator(&self, self.size());
    }
}

namespace Physica {
    template<class Derived>
    template<Vector V>
    class LValueVector<Derived>::View<V>::Iterator {
        using This = Iterator;
    public:
        using iterator_concept = std::contiguous_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = V::ScalarType;
        using reference = value_type::RefTy;
        using const_reference = value_type::ConstRefTy;
        using pointer = decltype(std::declval<V>().data_ptr(0));
    private:
        const LValueVector<Derived>::View<V>* view;
        difference_type index{};
    public:
        constexpr Iterator() = default;
        [[gnu::always_inline]] constexpr Iterator(const LValueVector<Derived>::View<V>* view, difference_type index) noexcept;
        constexpr Iterator(const This&) = default;
        constexpr Iterator(This&&) noexcept = default;
        constexpr ~Iterator() = default;
        /* Operators */
        constexpr This& operator=(const This&) = default;
        constexpr This& operator=(This&&) noexcept = default;
        [[gnu::always_inline]] constexpr This& operator++() noexcept;
        [[gnu::always_inline]] constexpr This& operator--() noexcept;
        [[gnu::always_inline]] constexpr This& operator+=(difference_type n) noexcept;
        [[gnu::always_inline]] constexpr This& operator-=(difference_type n) noexcept;
        [[nodiscard, gnu::always_inline]] constexpr This operator++(int) noexcept;
        [[nodiscard, gnu::always_inline]] constexpr This operator--(int) noexcept;
        [[nodiscard, gnu::always_inline]] constexpr decltype(auto) operator*() const noexcept;
        [[nodiscard, gnu::always_inline]] constexpr decltype(auto) operator[](difference_type n) const noexcept;
        [[nodiscard, gnu::always_inline]] constexpr bool operator==(const This& other) const noexcept;
        [[nodiscard, gnu::always_inline]] constexpr auto operator<=>(const This& other) const noexcept;
        [[nodiscard, gnu::always_inline]] constexpr This operator+(difference_type n) const noexcept;
        [[nodiscard, gnu::always_inline]] constexpr This operator-(difference_type n) const noexcept;
        [[nodiscard, gnu::always_inline]] constexpr difference_type operator-(const This& other) const noexcept;
        /* Operations */
        template<int Size>
        [[nodiscard, gnu::always_inline]] SIMD<value_type, Size> load() const noexcept;
        template<int Size>
        [[nodiscard, gnu::always_inline]] SIMD<value_type, Size> load(size_t count) const noexcept;
        [[gnu::always_inline]] void store(Packet auto pack) noexcept;
        [[gnu::always_inline]] void store(Packet auto pack, size_t count) noexcept;
        /* Friends */
        [[gnu::always_inline]] friend constexpr This operator+(difference_type n, const This& ite) noexcept { return ite + n; }
    };

    template<class Derived>
    template<Vector V>
    constexpr LValueVector<Derived>::View<V>::Iterator::Iterator(const LValueVector<Derived>::View<V>* view, difference_type index) noexcept : view(view), index(index) {}

    template<class Derived>
    template<Vector V>
    constexpr auto LValueVector<Derived>::View<V>::Iterator::operator++() noexcept -> This& {
        index += 1;
        return *this;
    }

    template<class Derived>
    template<Vector V>
    constexpr auto LValueVector<Derived>::View<V>::Iterator::operator--() noexcept -> This& {
        index -= 1;
        return *this;
    }

    template<class Derived>
    template<Vector V>
    constexpr auto LValueVector<Derived>::View<V>::Iterator::operator+=(difference_type n) noexcept -> This& {
        index += n;
        return *this;
    }

    template<class Derived>
    template<Vector V>
    constexpr auto LValueVector<Derived>::View<V>::Iterator::operator-=(difference_type n) noexcept -> This& {
        index -= n;
        return *this;
    }

    template<class Derived>
    template<Vector V>
    constexpr auto LValueVector<Derived>::View<V>::Iterator::operator++(int) noexcept -> This {
        return std::exchange(*this, This(view, index + 1));
    }

    template<class Derived>
    template<Vector V>
    constexpr auto LValueVector<Derived>::View<V>::Iterator::operator--(int) noexcept -> This {
        return std::exchange(*this, This(view, index - 1));
    }

    template<class Derived>
    template<Vector V>
    constexpr decltype(auto) LValueVector<Derived>::View<V>::Iterator::operator*() const noexcept {
        return view->vector()[index];
    }

    template<class Derived>
    template<Vector V>
    constexpr decltype(auto) LValueVector<Derived>::View<V>::Iterator::operator[](difference_type n) const noexcept {
        return view->vector()[index + n];
    }

    template<class Derived>
    template<Vector V>
    constexpr bool LValueVector<Derived>::View<V>::Iterator::operator==(const This& other) const noexcept {
        return index == other.index;
    }

    template<class Derived>
    template<Vector V>
    constexpr auto LValueVector<Derived>::View<V>::Iterator::operator<=>(const This& other) const noexcept {
        return index <=> other.index;
    }

    template<class Derived>
    template<Vector V>
    constexpr auto LValueVector<Derived>::View<V>::Iterator::operator+(difference_type n) const noexcept -> This {
        return This(view, index + n);
    }

    template<class Derived>
    template<Vector V>
    constexpr auto LValueVector<Derived>::View<V>::Iterator::operator-(difference_type n) const noexcept -> This {
        return This(view, index - n);
    }

    template<class Derived>
    template<Vector V>
    constexpr auto LValueVector<Derived>::View<V>::Iterator::operator-(const This& other) const noexcept -> difference_type {
        return index - other.index;
    }

    template<class Derived>
    template<Vector V>
    template<int Size>
    auto LValueVector<Derived>::View<V>::Iterator::load() const noexcept -> SIMD<value_type, Size> {
        return view->vector().template packet<Size>(index);
    }

    template<class Derived>
    template<Vector V>
    template<int Size>
    auto LValueVector<Derived>::View<V>::Iterator::load(size_t count) const noexcept -> SIMD<value_type, Size> {
        return view->vector().template packet<Size>(index, count);
    }

    template<class Derived>
    template<Vector V>
    void LValueVector<Derived>::View<V>::Iterator::store(Packet auto pack) noexcept {
        view->vector().writePacket(pack, index);
    }

    template<class Derived>
    template<Vector V>
    void LValueVector<Derived>::View<V>::Iterator::store(Packet auto pack, size_t count) noexcept {
        view->vector().writePacket(pack, index, count);
    }
}
