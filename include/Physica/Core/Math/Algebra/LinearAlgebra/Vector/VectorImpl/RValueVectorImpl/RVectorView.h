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

#include "../RValueVector.h"

namespace Physica {
    /**
     * Vector views are intended for loop-invariant code motion. We do this because compilers struggle to prove noalias and to align with the sprit of C++20 ranges.
     */
    template<Vector V>
    class RVectorView : public std::ranges::view_base {
        static_assert(!std::is_reference_v<V>);
        using This = RVectorView<V>;

        class Iterator;

        V* vec;
    public:
        [[gnu::always_inline]] constexpr RVectorView(V& vec) noexcept : vec(&vec) {}
        constexpr RVectorView(const This&) = default;
        constexpr RVectorView(This&&) noexcept = default;
        constexpr ~RVectorView() = default;
        /* Operators */
        constexpr This& operator=(const This&) = default;
        constexpr This& operator=(This&&) noexcept = default;
        /* Operations */
        [[nodiscard, gnu::always_inline]] decltype(auto) calc(size_t index) const noexcept;

        [[nodiscard, gnu::always_inline]] constexpr auto begin(this auto&&) noexcept;
        [[nodiscard, gnu::always_inline]] constexpr auto end(this auto&&) noexcept;
        /* Getters */
        [[nodiscard, gnu::always_inline]] constexpr V& vector() const noexcept { return *vec; }
        [[nodiscard, gnu::always_inline]] constexpr size_t size() const noexcept;
    };

    template<Vector V>
    decltype(auto) RVectorView<V>::calc(size_t index) const noexcept {
        return vec->calc(index);
    }

    template<Vector V>
    constexpr auto RVectorView<V>::begin(this auto&& self) noexcept {
        return Iterator(&self, 0);
    }

    template<Vector V>
    constexpr auto RVectorView<V>::end(this auto&& self) noexcept {
        return Iterator(&self, self.size());
    }

    template<Vector V>
    constexpr size_t RVectorView<V>::size() const noexcept {
        return vector().getLength();
    }
}

namespace Physica {
    template<Vector V>
    class RVectorView<V>::Iterator {
        using This = Iterator;
    public:
        using iterator_concept = std::random_access_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = V::ScalarType;
        using reference = const value_type;
        using const_reference = const value_type;
    private:
        const RVectorView<V>* view;
        difference_type index{};
    public:
        constexpr Iterator() = default;
        [[gnu::always_inline]] constexpr Iterator(const RVectorView<V>* view, difference_type index) noexcept;
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
        [[nodiscard, gnu::always_inline]] constexpr reference operator*() const noexcept;
        [[nodiscard, gnu::always_inline]] constexpr reference operator[](difference_type n) const noexcept;
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
        /* Friends */
        [[gnu::always_inline]] friend constexpr This operator+(difference_type n, const This& ite) noexcept { return ite + n; }
    };

    template<Vector V>
    constexpr RVectorView<V>::Iterator::Iterator(const RVectorView<V>* view, difference_type index) noexcept : view(view), index(index) {}

    template<Vector V>
    constexpr auto RVectorView<V>::Iterator::operator++() noexcept -> This& {
        index += 1;
        return *this;
    }

    template<Vector V>
    constexpr auto RVectorView<V>::Iterator::operator--() noexcept -> This& {
        index -= 1;
        return *this;
    }

    template<Vector V>
    constexpr auto RVectorView<V>::Iterator::operator+=(difference_type n) noexcept -> This& {
        index += n;
        return *this;
    }

    template<Vector V>
    constexpr auto RVectorView<V>::Iterator::operator-=(difference_type n) noexcept -> This& {
        index -= n;
        return *this;
    }

    template<Vector V>
    constexpr auto RVectorView<V>::Iterator::operator++(int) noexcept -> This {
        return std::exchange(*this, This(view, index + 1));
    }

    template<Vector V>
    constexpr auto RVectorView<V>::Iterator::operator--(int) noexcept -> This {
        return std::exchange(*this, This(view, index - 1));
    }

    template<Vector V>
    constexpr auto RVectorView<V>::Iterator::operator*() const noexcept -> reference {
        return view->calc(index);
    }

    template<Vector V>
    constexpr auto RVectorView<V>::Iterator::operator[](difference_type n) const noexcept -> reference {
        return view->calc(index + n);
    }

    template<Vector V>
    constexpr bool RVectorView<V>::Iterator::operator==(const This& other) const noexcept {
        return index == other.index;
    }

    template<Vector V>
    constexpr auto RVectorView<V>::Iterator::operator<=>(const This& other) const noexcept {
        return index <=> other.index;
    }

    template<Vector V>
    constexpr auto RVectorView<V>::Iterator::operator+(difference_type n) const noexcept -> This {
        return This(view, index + n);
    }

    template<Vector V>
    constexpr auto RVectorView<V>::Iterator::operator-(difference_type n) const noexcept -> This {
        return This(view, index - n);
    }

    template<Vector V>
    constexpr auto RVectorView<V>::Iterator::operator-(const This& other) const noexcept -> difference_type {
        return index - other.index;
    }

    template<Vector V>
    template<int Size>
    auto RVectorView<V>::Iterator::load() const noexcept -> SIMD<value_type, Size> {
        assert(index + Size <= view->size());
        return view->vector().template packet<Size>(index);
    }

    template<Vector V>
    template<int Size>
    auto RVectorView<V>::Iterator::load(size_t count) const noexcept -> SIMD<value_type, Size> {
        assert(0 < count && count < Size && "[Error]: Invalid size for partial operation");
        return view->vector().template packet<Size>(index, count);
    }
}
