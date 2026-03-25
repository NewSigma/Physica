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
    template<ExprID ID, class LHS, class RHS>
    template<class ViewLHS, class ViewRHS>
    class BinaryVectorExpr<ID, LHS, RHS>::View : public std::ranges::view_base {
        using This = View<ViewLHS, ViewRHS>;

        class Iterator;

        ViewLHS lhs;
        ViewRHS rhs;
    public:
        [[gnu::always_inline]] constexpr View(ViewLHS lhs, ViewRHS rhs);
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
        [[nodiscard, gnu::always_inline]] constexpr size_t size() const noexcept { return lhs.size(); }
    };

    template<ExprID ID, class LHS, class RHS>
    template<class ViewLHS, class ViewRHS>
    constexpr BinaryVectorExpr<ID, LHS, RHS>::View<ViewLHS, ViewRHS>::View(ViewLHS lhs, ViewRHS rhs)
            : lhs(std::move(lhs)), rhs(std::move(rhs)) {}

    template<ExprID ID, class LHS, class RHS>
    template<class ViewLHS, class ViewRHS>
    constexpr auto BinaryVectorExpr<ID, LHS, RHS>::View<ViewLHS, ViewRHS>::begin(this auto&& self) noexcept {
        if constexpr (Scalar<LHS>)
            return Iterator(self.lhs, self.rhs.begin());
        else if constexpr (Scalar<RHS>)
            return Iterator(self.lhs.begin(), self.rhs);
        else
            return Iterator(self.lhs.begin(), self.rhs.begin());
    }

    template<ExprID ID, class LHS, class RHS>
    template<class ViewLHS, class ViewRHS>
    constexpr auto BinaryVectorExpr<ID, LHS, RHS>::View<ViewLHS, ViewRHS>::end(this auto&& self) noexcept {
        if constexpr (Scalar<LHS>)
            return Iterator(self.lhs, self.rhs.end());
        else if constexpr (Scalar<RHS>)
            return Iterator(self.lhs.end(), self.rhs);
        else
            return Iterator(self.lhs.end(), self.rhs.end());
    }
}

namespace Physica {
    template<ExprID ID, class LHS, class RHS>
    template<class ViewLHS, class ViewRHS>
    class BinaryVectorExpr<ID, LHS, RHS>::View<ViewLHS, ViewRHS>::Iterator {
        using This = Iterator;
        template<class T>
        struct IteratorOrScalar {
            using Type = std::ranges::iterator_t<T>;
        };

        template<Scalar T>
        struct IteratorOrScalar<T> {
            using Type = T;
        };

        using ItLHS = IteratorOrScalar<ViewLHS>::Type;
        using ItRHS = IteratorOrScalar<ViewRHS>::Type;
    public:
        using iterator_concept = std::contiguous_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = VectorExpr<ID, LHS, RHS>::ScalarType;
        using reference = const value_type;
        using const_reference = const value_type;
    private:
        ItLHS lhs;
        ItRHS rhs;
    public:
        constexpr Iterator() = default;
        [[gnu::always_inline]] constexpr Iterator(ItLHS lhs, ItRHS rhs) noexcept;
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
        /* Getters */
        [[nodiscard, gnu::always_inline]] auto getLHS() const noexcept { return lhs; }
        [[nodiscard, gnu::always_inline]] auto getRHS() const noexcept { return rhs; }
        /* Friends */
        [[gnu::always_inline]] friend constexpr This operator+(difference_type n, const This& ite) noexcept { return ite + n; }
    };

    template<ExprID ID, class LHS, class RHS>
    template<class ViewLHS, class ViewRHS>
    constexpr BinaryVectorExpr<ID, LHS, RHS>::View<ViewLHS, ViewRHS>::Iterator::Iterator(ItLHS lhs, ItRHS rhs) noexcept : lhs(lhs), rhs(rhs) {}

    template<ExprID ID, class LHS, class RHS>
    template<class ViewLHS, class ViewRHS>
    constexpr auto BinaryVectorExpr<ID, LHS, RHS>::View<ViewLHS, ViewRHS>::Iterator::operator++() noexcept -> This& {
        if constexpr (!Scalar<ItLHS>)
            ++lhs;
        if constexpr (!Scalar<ItRHS>)
            ++rhs;
        return *this;
    }

    template<ExprID ID, class LHS, class RHS>
    template<class ViewLHS, class ViewRHS>
    constexpr auto BinaryVectorExpr<ID, LHS, RHS>::View<ViewLHS, ViewRHS>::Iterator::operator--() noexcept -> This& {
        if constexpr (!Scalar<ItLHS>)
            --lhs;
        if constexpr (!Scalar<ItRHS>)
            --rhs;
        return *this;
    }

    template<ExprID ID, class LHS, class RHS>
    template<class ViewLHS, class ViewRHS>
    constexpr auto BinaryVectorExpr<ID, LHS, RHS>::View<ViewLHS, ViewRHS>::Iterator::operator+=(difference_type n) noexcept -> This& {
        if constexpr (!Scalar<ItLHS>)
            lhs += n;
        if constexpr (!Scalar<ItRHS>)
            rhs += n;
        return *this;
    }

    template<ExprID ID, class LHS, class RHS>
    template<class ViewLHS, class ViewRHS>
    constexpr auto BinaryVectorExpr<ID, LHS, RHS>::View<ViewLHS, ViewRHS>::Iterator::operator-=(difference_type n) noexcept -> This& {
        if constexpr (!Scalar<ItLHS>)
            lhs -= n;
        if constexpr (!Scalar<ItRHS>)
            rhs -= n;
        return *this;
    }

    template<ExprID ID, class LHS, class RHS>
    template<class ViewLHS, class ViewRHS>
    constexpr auto BinaryVectorExpr<ID, LHS, RHS>::View<ViewLHS, ViewRHS>::Iterator::operator++(int) noexcept -> This {
        This temp = *this;
        if constexpr (!Scalar<ItLHS>)
            ++lhs;
        if constexpr (!Scalar<ItRHS>)
            ++rhs;
        return temp;
    }

    template<ExprID ID, class LHS, class RHS>
    template<class ViewLHS, class ViewRHS>
    constexpr auto BinaryVectorExpr<ID, LHS, RHS>::View<ViewLHS, ViewRHS>::Iterator::operator--(int) noexcept -> This {
        This temp = *this;
        if constexpr (!Scalar<ItLHS>)
            --lhs;
        if constexpr (!Scalar<ItRHS>)
            --rhs;
        return temp;
    }

    template<ExprID ID, class LHS, class RHS>
    template<class ViewLHS, class ViewRHS>
    constexpr auto BinaryVectorExpr<ID, LHS, RHS>::View<ViewLHS, ViewRHS>::Iterator::operator*() const noexcept -> reference {
        return VectorExpr<ID, LHS, RHS>::operator()(lhs, rhs);
    }

    template<ExprID ID, class LHS, class RHS>
    template<class ViewLHS, class ViewRHS>
    constexpr auto BinaryVectorExpr<ID, LHS, RHS>::View<ViewLHS, ViewRHS>::Iterator::operator[](difference_type n) const noexcept -> reference {
        if constexpr (Scalar<ItLHS>)
            return VectorExpr<ID, LHS, RHS>::operator()(lhs, rhs + n);
        else if constexpr (Scalar<ItRHS>)
            return VectorExpr<ID, LHS, RHS>::operator()(lhs + n, rhs);
        else
            return VectorExpr<ID, LHS, RHS>::operator()(lhs + n, rhs + n);
    }

    template<ExprID ID, class LHS, class RHS>
    template<class ViewLHS, class ViewRHS>
    constexpr bool BinaryVectorExpr<ID, LHS, RHS>::View<ViewLHS, ViewRHS>::Iterator::operator==(const This& other) const noexcept {
        if constexpr (Scalar<ItLHS>)
            return rhs == other.rhs;
        else
            return lhs == other.lhs;
    }

    template<ExprID ID, class LHS, class RHS>
    template<class ViewLHS, class ViewRHS>
    constexpr auto BinaryVectorExpr<ID, LHS, RHS>::View<ViewLHS, ViewRHS>::Iterator::operator<=>(const This& other) const noexcept {
        if constexpr (Scalar<ItLHS>)
            return rhs <=> other.rhs;
        else
            return lhs <=> other.lhs;
    }

    template<ExprID ID, class LHS, class RHS>
    template<class ViewLHS, class ViewRHS>
    constexpr auto BinaryVectorExpr<ID, LHS, RHS>::View<ViewLHS, ViewRHS>::Iterator::operator+(difference_type n) const noexcept -> This {
        if constexpr (Scalar<ItLHS>)
            return This(lhs, rhs + n);
        else if constexpr (Scalar<ItRHS>)
            return This(lhs + n, rhs);
        else
            return This(lhs + n, rhs + n);
    }

    template<ExprID ID, class LHS, class RHS>
    template<class ViewLHS, class ViewRHS>
    constexpr auto BinaryVectorExpr<ID, LHS, RHS>::View<ViewLHS, ViewRHS>::Iterator::operator-(difference_type n) const noexcept -> This {
        if constexpr (Scalar<ItLHS>)
            return This(lhs, rhs - n);
        else if constexpr (Scalar<ItRHS>)
            return This(lhs - n, rhs);
        else
            return This(lhs - n, rhs - n);
    }

    template<ExprID ID, class LHS, class RHS>
    template<class ViewLHS, class ViewRHS>
    constexpr auto BinaryVectorExpr<ID, LHS, RHS>::View<ViewLHS, ViewRHS>::Iterator::operator-(const This& other) const noexcept -> difference_type {
        if constexpr (Scalar<ItLHS>)
            return rhs - other.rhs;
        else
            return lhs - other.lhs;
    }

    template<ExprID ID, class LHS, class RHS>
    template<class ViewLHS, class ViewRHS>
    template<int Size>
    auto BinaryVectorExpr<ID, LHS, RHS>::View<ViewLHS, ViewRHS>::Iterator::load() const noexcept -> SIMD<value_type, Size> {
        return VectorExpr<ID, LHS, RHS>::template operator()<Size>(lhs, rhs);
    }

    template<ExprID ID, class LHS, class RHS>
    template<class ViewLHS, class ViewRHS>
    template<int Size>
    auto BinaryVectorExpr<ID, LHS, RHS>::View<ViewLHS, ViewRHS>::Iterator::load(size_t count) const noexcept -> SIMD<value_type, Size> {
        return VectorExpr<ID, LHS, RHS>::template operator()<Size>(lhs, rhs, count);
    }
}
