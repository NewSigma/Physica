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

#include "../CompactVector.h"

namespace Physica {
    template<class Derived>
    template<Vector V>
    class CompactVector<Derived>::View final : public std::ranges::view_base {
        using This = View;

        using pointer = decltype(std::declval<V>().data());

        class Iterator;
    private:
        pointer data;
        size_t length;
    public:
        [[gnu::always_inline]] constexpr View(V& vec) noexcept;
        constexpr View(const This&) = default;
        constexpr View(This&&) noexcept = default;
        constexpr ~View() = default;
        /* Operators */
        constexpr This& operator=(const This&) = default;
        constexpr This& operator=(This&&) noexcept = default;
        /* Operations */
        [[nodiscard, gnu::always_inline]] constexpr auto begin(this auto&&) noexcept;
        [[nodiscard, gnu::always_inline]] constexpr auto end(this auto&&) noexcept;
        [[nodiscard, gnu::always_inline]] constexpr size_t size() const noexcept;
    };

    template<class Derived>
    template<Vector V>
    constexpr CompactVector<Derived>::View<V>::View(V& vec) noexcept : data(vec.data()), length(vec.getLength()) {}

    template<class Derived>
    template<Vector V>
    constexpr auto CompactVector<Derived>::View<V>::begin(this auto&& self) noexcept {
        return Iterator(self.data);
    }

    template<class Derived>
    template<Vector V>
    constexpr auto CompactVector<Derived>::View<V>::end(this auto&& self) noexcept {
        return Iterator(self.data + self.length);
    }

    template<class Derived>
    template<Vector V>
    constexpr size_t CompactVector<Derived>::View<V>::size() const noexcept {
        return length;
    }
}

namespace Physica {
    template<class Derived>
    template<Vector V>
    class CompactVector<Derived>::View<V>::Iterator {
        using This = Iterator;
    public:
        using iterator_concept = std::contiguous_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = V::ScalarType;
        using reference = value_type::RefTy;
        using const_reference = value_type::ConstRefTy;
        using pointer = decltype(std::declval<V>().data());
    private:
        pointer pos;
    public:
        constexpr Iterator() = default;
        [[gnu::always_inline]] constexpr Iterator(pointer pos) noexcept;
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
    constexpr CompactVector<Derived>::View<V>::Iterator::Iterator(pointer pos) noexcept : pos(pos) {}

    template<class Derived>
    template<Vector V>
    constexpr auto CompactVector<Derived>::View<V>::Iterator::operator++() noexcept -> This& {
        pos += 1;
        return *this;
    }

    template<class Derived>
    template<Vector V>
    constexpr auto CompactVector<Derived>::View<V>::Iterator::operator--() noexcept -> This& {
        pos -= 1;
        return *this;
    }

    template<class Derived>
    template<Vector V>
    constexpr auto CompactVector<Derived>::View<V>::Iterator::operator+=(difference_type n) noexcept -> This& {
        pos += n;
        return *this;
    }

    template<class Derived>
    template<Vector V>
    constexpr auto CompactVector<Derived>::View<V>::Iterator::operator-=(difference_type n) noexcept -> This& {
        pos -= n;
        return *this;
    }

    template<class Derived>
    template<Vector V>
    constexpr auto CompactVector<Derived>::View<V>::Iterator::operator++(int) noexcept -> This {
        return std::exchange(*this, This(pos + 1));
    }

    template<class Derived>
    template<Vector V>
    constexpr auto CompactVector<Derived>::View<V>::Iterator::operator--(int) noexcept -> This {
        return std::exchange(*this, This(pos - 1));
    }

    template<class Derived>
    template<Vector V>
    constexpr decltype(auto) CompactVector<Derived>::View<V>::Iterator::operator*() const noexcept {
        return *pos;
    }

    template<class Derived>
    template<Vector V>
    constexpr decltype(auto) CompactVector<Derived>::View<V>::Iterator::operator[](difference_type n) const noexcept {
        return pos[n];
    }

    template<class Derived>
    template<Vector V>
    constexpr bool CompactVector<Derived>::View<V>::Iterator::operator==(const This& other) const noexcept {
        return pos == other.pos;
    }

    template<class Derived>
    template<Vector V>
    constexpr auto CompactVector<Derived>::View<V>::Iterator::operator<=>(const This& other) const noexcept {
        return pos <=> other.pos;
    }

    template<class Derived>
    template<Vector V>
    constexpr auto CompactVector<Derived>::View<V>::Iterator::operator+(difference_type n) const noexcept -> This {
        return This(pos + n);
    }

    template<class Derived>
    template<Vector V>
    constexpr auto CompactVector<Derived>::View<V>::Iterator::operator-(difference_type n) const noexcept -> This {
        return This(pos - n);
    }

    template<class Derived>
    template<Vector V>
    constexpr auto CompactVector<Derived>::View<V>::Iterator::operator-(const This& other) const noexcept -> difference_type {
        return pos - other.pos;
    }

    template<class Derived>
    template<Vector V>
    template<int Size>
    auto CompactVector<Derived>::View<V>::Iterator::load() const noexcept -> SIMD<value_type, Size> {
        SIMD<value_type, Size> pack{};
        pack.load(pos);
        return pack;
    }

    template<class Derived>
    template<Vector V>
    template<int Size>
    auto CompactVector<Derived>::View<V>::Iterator::load(size_t count) const noexcept -> SIMD<value_type, Size> {
        assert(0 < count && count < Size && "[Error]: Invalid size for partial operation");
        SIMD<value_type, Size> pack{};
        pack.load(pos, count);
        return pack;
    }

    template<class Derived>
    template<Vector V>
    void CompactVector<Derived>::View<V>::Iterator::store(const Packet auto pack) noexcept {
        pack.store(pos);
    }

    template<class Derived>
    template<Vector V>
    void CompactVector<Derived>::View<V>::Iterator::store(const Packet auto pack, size_t count) noexcept {
        assert(0 < count && count < pack.size() && "[Error]: Invalid size for partial operation");
        pack.store(pos, count);
    }
}
