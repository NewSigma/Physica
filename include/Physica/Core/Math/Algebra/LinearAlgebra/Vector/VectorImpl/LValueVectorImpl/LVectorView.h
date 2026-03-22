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
    template<Vector V>
    class LVectorView : public RVectorView<V> {
        using This = LVectorView<V>;
        using Base = RVectorView<V>;

        class Iterator;
    public:
        constexpr LVectorView(V& vec) noexcept : Base(vec) {}
        constexpr LVectorView(const This&) = default;
        constexpr LVectorView(This&&) noexcept = default;
        constexpr ~LVectorView() = default;
        /* Operators */
        constexpr This& operator=(const This&) = default;
        constexpr This& operator=(This&&) noexcept = default;
        /* Operations */
        [[nodiscard]] constexpr auto begin(this auto&&) noexcept;
        [[nodiscard]] constexpr auto end(this auto&&) noexcept;
    };

    template<Vector V>
    constexpr auto LVectorView<V>::begin(this auto&& self) noexcept {
        return Iterator(&self, 0);
    }

    template<Vector V>
    constexpr auto LVectorView<V>::end(this auto&& self) noexcept {
        return Iterator(&self, self.size());
    }
}

namespace Physica {
    template<Vector V>
    class LVectorView<V>::Iterator {
        using This = Iterator;
    public:
        using iterator_concept = std::contiguous_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = V::ScalarType;
        using reference = value_type::RefTy;
        using const_reference = value_type::ConstRefTy;
        using pointer = std::conditional<std::is_const_v<V>, typename value_type::ConstPtrTy, typename value_type::PtrTy>::type;
    private:
        const LVectorView<V>* view;
        difference_type index{};
    public:
        constexpr Iterator() = default;
        constexpr Iterator(const LVectorView<V>* view, difference_type index) noexcept;
        constexpr Iterator(const This&) = default;
        constexpr Iterator(This&&) noexcept = default;
        constexpr ~Iterator() = default;
        /* Operators */
        constexpr This& operator=(const This&) = default;
        constexpr This& operator=(This&&) noexcept = default;
        constexpr This& operator++() noexcept;
        constexpr This& operator--() noexcept;
        constexpr This& operator+=(difference_type n) noexcept;
        constexpr This& operator-=(difference_type n) noexcept;
        [[nodiscard]] constexpr This operator++(int) noexcept;
        [[nodiscard]] constexpr This operator--(int) noexcept;
        [[nodiscard]] constexpr decltype(auto) operator*() const noexcept;
        [[nodiscard]] constexpr decltype(auto) operator[](difference_type n) const noexcept;
        [[nodiscard]] constexpr bool operator==(const This& other) const noexcept;
        [[nodiscard]] constexpr auto operator<=>(const This& other) const noexcept;
        [[nodiscard]] constexpr This operator+(difference_type n) const noexcept;
        [[nodiscard]] constexpr This operator-(difference_type n) const noexcept;
        [[nodiscard]] constexpr difference_type operator-(const This& other) const noexcept;
        /* Operations */
        template<int Size>
        [[nodiscard]] SIMD<value_type, Size> load() const noexcept;
        template<int Size>
        [[nodiscard]] SIMD<value_type, Size> load(size_t count) const noexcept;
        void store(Packet auto pack) noexcept;
        void store(Packet auto pack, size_t count) noexcept;
        /* Friends */
        friend constexpr This operator+(difference_type n, const This& ite) noexcept { return ite + n; }
    };

    template<Vector V>
    constexpr LVectorView<V>::Iterator::Iterator(const LVectorView<V>* view, difference_type index) noexcept : view(view), index(index) {}

    template<Vector V>
    constexpr auto LVectorView<V>::Iterator::operator++() noexcept -> This& {
        index += 1;
        return *this;
    }

    template<Vector V>
    constexpr auto LVectorView<V>::Iterator::operator--() noexcept -> This& {
        index -= 1;
        return *this;
    }

    template<Vector V>
    constexpr auto LVectorView<V>::Iterator::operator+=(difference_type n) noexcept -> This& {
        index += n;
        return *this;
    }

    template<Vector V>
    constexpr auto LVectorView<V>::Iterator::operator-=(difference_type n) noexcept -> This& {
        index -= n;
        return *this;
    }

    template<Vector V>
    constexpr auto LVectorView<V>::Iterator::operator++(int) noexcept -> This {
        return std::exchange(*this, This(view, index + 1));
    }

    template<Vector V>
    constexpr auto LVectorView<V>::Iterator::operator--(int) noexcept -> This {
        return std::exchange(*this, This(view, index - 1));
    }

    template<Vector V>
    constexpr decltype(auto) LVectorView<V>::Iterator::operator*() const noexcept {
        return view->vector()[index];
    }

    template<Vector V>
    constexpr decltype(auto) LVectorView<V>::Iterator::operator[](difference_type n) const noexcept {
        return view->vector()[index + n];
    }

    template<Vector V>
    constexpr bool LVectorView<V>::Iterator::operator==(const This& other) const noexcept {
        return index == other.index;
    }

    template<Vector V>
    constexpr auto LVectorView<V>::Iterator::operator<=>(const This& other) const noexcept {
        return index <=> other.index;
    }

    template<Vector V>
    constexpr auto LVectorView<V>::Iterator::operator+(difference_type n) const noexcept -> This {
        return This(view, index + n);
    }

    template<Vector V>
    constexpr auto LVectorView<V>::Iterator::operator-(difference_type n) const noexcept -> This {
        return This(view, index - n);
    }

    template<Vector V>
    constexpr auto LVectorView<V>::Iterator::operator-(const This& other) const noexcept -> difference_type {
        return index - other.index;
    }

    template<Vector V>
    template<int Size>
    auto LVectorView<V>::Iterator::load() const noexcept -> SIMD<value_type, Size> {
        return view->vector().template packet<Size>(index);
    }

    template<Vector V>
    template<int Size>
    auto LVectorView<V>::Iterator::load(size_t count) const noexcept -> SIMD<value_type, Size> {
        return view->vector().template packet<Size>(index, count);
    }

    template<Vector V>
    void LVectorView<V>::Iterator::store(Packet auto pack) noexcept {
        view->vector().writePacket(pack, index);
    }

    template<Vector V>
    void LVectorView<V>::Iterator::store(Packet auto pack, size_t count) noexcept {
        view->vector().writePacket(pack, index, count);
    }
}
