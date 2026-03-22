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

#include "../VectorExpr.h"

namespace Physica {
    template<ExprID ID, Vector V>
    template<std::ranges::view Operand>
    class UnitaryVectorExpr<ID, V>::View : public std::ranges::view_base {
        using This = View<Operand>;

        class Iterator;

        Operand operand;
    public:
        constexpr View(Operand operand);
        constexpr View(const This&) = default;
        constexpr View(This&&) noexcept = default;
        constexpr ~View() = default;
        /* Operators */
        constexpr This& operator=(const This&) = default;
        constexpr This& operator=(This&&) noexcept = default;
        /* Operations */
        [[nodiscard]] constexpr auto begin(this auto&&) noexcept;
        [[nodiscard]] constexpr auto end(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] constexpr size_t size() const noexcept { return operand.size(); }
    };

    template<ExprID ID, Vector V>
    template<std::ranges::view Operand>
    constexpr UnitaryVectorExpr<ID, V>::View<Operand>::View(Operand operand) : operand(std::move(operand)) {}

    template<ExprID ID, Vector V>
    template<std::ranges::view Operand>
    constexpr auto UnitaryVectorExpr<ID, V>::View<Operand>::begin(this auto&& self) noexcept {
        return Iterator(self.operand.begin());
    }

    template<ExprID ID, Vector V>
    template<std::ranges::view Operand>
    constexpr auto UnitaryVectorExpr<ID, V>::View<Operand>::end(this auto&& self) noexcept {
        return Iterator(self.operand.end());
    }
}

namespace Physica {
    template<ExprID ID, Vector V>
    template<std::ranges::view Operand>
    class UnitaryVectorExpr<ID, V>::View<Operand>::Iterator {
        using This = Iterator;
    public:
        using iterator_concept = std::random_access_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = VectorExpr<ID, V>::ScalarType;
        using reference = const value_type;
        using const_reference = const value_type;
    private:
        std::ranges::iterator_t<Operand> it;
    public:
        constexpr Iterator() = default;
        constexpr Iterator(std::ranges::iterator_t<Operand> it) noexcept;
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
        [[nodiscard]] constexpr reference operator*() const noexcept;
        [[nodiscard]] constexpr reference operator[](difference_type n) const noexcept;
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
        /* Friends */
        friend constexpr This operator+(difference_type n, const This& ite) noexcept { return ite + n; }
    };

    template<ExprID ID, Vector V>
    template<std::ranges::view Operand>
    constexpr UnitaryVectorExpr<ID, V>::View<Operand>::Iterator::Iterator(std::ranges::iterator_t<Operand> it) noexcept : it(std::move(it)) {}

    template<ExprID ID, Vector V>
    template<std::ranges::view Operand>
    constexpr auto UnitaryVectorExpr<ID, V>::View<Operand>::Iterator::operator++() noexcept -> This& {
        ++it;
        return *this;
    }

    template<ExprID ID, Vector V>
    template<std::ranges::view Operand>
    constexpr auto UnitaryVectorExpr<ID, V>::View<Operand>::Iterator::operator--() noexcept -> This& {
        --it;
        return *this;
    }

    template<ExprID ID, Vector V>
    template<std::ranges::view Operand>
    constexpr auto UnitaryVectorExpr<ID, V>::View<Operand>::Iterator::operator+=(difference_type n) noexcept -> This& {
        it += n;
        return *this;
    }

    template<ExprID ID, Vector V>
    template<std::ranges::view Operand>
    constexpr auto UnitaryVectorExpr<ID, V>::View<Operand>::Iterator::operator-=(difference_type n) noexcept -> This& {
        it -= n;
        return *this;
    }

    template<ExprID ID, Vector V>
    template<std::ranges::view Operand>
    constexpr auto UnitaryVectorExpr<ID, V>::View<Operand>::Iterator::operator++(int) noexcept -> This {
        return std::exchange(*this, it++);
    }

    template<ExprID ID, Vector V>
    template<std::ranges::view Operand>
    constexpr auto UnitaryVectorExpr<ID, V>::View<Operand>::Iterator::operator--(int) noexcept -> This {
        return std::exchange(*this, it--);
    }

    template<ExprID ID, Vector V>
    template<std::ranges::view Operand>
    constexpr auto UnitaryVectorExpr<ID, V>::View<Operand>::Iterator::operator*() const noexcept -> reference {
        return VectorExpr<ID, V>::operator()(it);
    }

    template<ExprID ID, Vector V>
    template<std::ranges::view Operand>
    constexpr auto UnitaryVectorExpr<ID, V>::View<Operand>::Iterator::operator[](difference_type n) const noexcept -> reference {
        return VectorExpr<ID, V>::operator()(it + n);
    }

    template<ExprID ID, Vector V>
    template<std::ranges::view Operand>
    constexpr bool UnitaryVectorExpr<ID, V>::View<Operand>::Iterator::operator==(const This& other) const noexcept {
        return it == other.it;
    }

    template<ExprID ID, Vector V>
    template<std::ranges::view Operand>
    constexpr auto UnitaryVectorExpr<ID, V>::View<Operand>::Iterator::operator<=>(const This& other) const noexcept {
        return it <=> other.it;
    }

    template<ExprID ID, Vector V>
    template<std::ranges::view Operand>
    constexpr auto UnitaryVectorExpr<ID, V>::View<Operand>::Iterator::operator+(difference_type n) const noexcept -> This {
        return This(it + n);
    }

    template<ExprID ID, Vector V>
    template<std::ranges::view Operand>
    constexpr auto UnitaryVectorExpr<ID, V>::View<Operand>::Iterator::operator-(difference_type n) const noexcept -> This {
        return This(it - n);
    }

    template<ExprID ID, Vector V>
    template<std::ranges::view Operand>
    constexpr auto UnitaryVectorExpr<ID, V>::View<Operand>::Iterator::operator-(const This& other) const noexcept -> difference_type {
        return it - other.it;
    }

    template<ExprID ID, Vector V>
    template<std::ranges::view Operand>
    template<int Size>
    auto UnitaryVectorExpr<ID, V>::View<Operand>::Iterator::load() const noexcept -> SIMD<value_type, Size> {
        return VectorExpr<ID, V>::template operator()<Size>(it);
    }

    template<ExprID ID, Vector V>
    template<std::ranges::view Operand>
    template<int Size>
    auto UnitaryVectorExpr<ID, V>::View<Operand>::Iterator::load(size_t count) const noexcept -> SIMD<value_type, Size> {
        return VectorExpr<ID, V>::template operator()<Size>(it, count);
    }
}
