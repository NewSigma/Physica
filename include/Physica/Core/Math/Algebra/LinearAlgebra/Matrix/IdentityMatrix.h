/*
 * Copyright 2024-2026 Weibo He.
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

#include "MatrixImpl/RValueMatrix.h"

namespace Physica {
    template<Scalar T, size_t Order = Dynamic>
    class IdentityMatrix : public RValueMatrix<IdentityMatrix<T, Order>> {
        using This = IdentityMatrix<T, Order>;
        using Base = RValueMatrix<This>;
        using IndexType = std::conditional<Order == Dynamic, size_t, PlainStruct<void>>::type;
        static_assert(!T::isComplex && !T::isDiffable, "[Error]: Invalid scalar for unit matrix");
    protected:
        using typename Base::Tv;
    private:
        [[no_unique_address]] IndexType order;
    public:
        IdentityMatrix() = default;
        explicit IdentityMatrix(size_t order_);
        explicit IdentityMatrix(const Matrix auto& m);
        IdentityMatrix(const This&) = default;
        IdentityMatrix(This&&) noexcept = default;
        ~IdentityMatrix() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] auto&& operator*(Vector auto&& v) const noexcept;
        [[nodiscard]] auto&& operator*(Matrix auto&& m) const noexcept;
        /* Operations */
        void assign(Matrix auto&& target) const;

        [[nodiscard]] T calc(size_t row, size_t col) const;
        [[nodiscard]] T calc_value(size_t row, size_t col) const;

        [[nodiscard]] const This& transpose() const noexcept { return *this; }
        [[nodiscard]] const This& conjugate() const noexcept { return *this; }
        [[nodiscard]] const This& hermite() const noexcept { return *this; }
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept;
        [[nodiscard]] size_t getCol() const noexcept;
        /* Friends */
        friend class device_obj<This>;
    };

    template<Scalar T, size_t Order>
    IdentityMatrix<T, Order>::IdentityMatrix(size_t order_) : order(order_) {
        assert((Order == Dynamic || Order == order_) && "[Error]: tparam and param is not consistent");
    }

    template<Scalar T, size_t Order>
    IdentityMatrix<T, Order>::IdentityMatrix(const Matrix auto& m) : IdentityMatrix(m.getRow()) {
        assert(m.isSquare());
    }

    template<Scalar T, size_t Order>
    auto&& IdentityMatrix<T, Order>::operator*(Vector auto&& v) const noexcept {
        using RHS = decltype(v);
        assert(getCol() == v.getLength());
        return std::forward<RHS>(v);
    }

    template<Scalar T, size_t Order>
    auto&& IdentityMatrix<T, Order>::operator*(Matrix auto&& m) const noexcept {
        using RHS = decltype(m);
        assert(getCol() == m.getRow());
        return std::forward<RHS>(m);
    }

    template<Scalar T, size_t Order>
    void IdentityMatrix<T, Order>::assign(Matrix auto&& target) const {
        target.zeros();
        for (size_t i = 0; i < getRow(); ++i)
            target(i, i) = T(1);
    }

    template<Scalar T, size_t Order>
    T IdentityMatrix<T, Order>::calc(size_t row, size_t col) const {
        return T(row == col ? 1 : 0);
    }

    template<Scalar T, size_t Order>
    T IdentityMatrix<T, Order>::calc_value(size_t row, size_t col) const {
        return calc(row, col);
    }

    template<Scalar T, size_t Order>
    void IdentityMatrix<T, Order>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(order, obj.order);
    }

    template<Scalar T, size_t Order>
    size_t IdentityMatrix<T, Order>::getRow() const noexcept {
        if constexpr (Order == Dynamic)
            return order;
        return Order;
    }

    template<Scalar T, size_t Order>
    size_t IdentityMatrix<T, Order>::getCol() const noexcept {
        return getRow();
    }
}

namespace Physica {
    template<Scalar T, size_t Order>
    class Traits<IdentityMatrix<T, Order>> {
    public:
        using ScalarType = T;
        constexpr static int Major = MatrixMajor::BothMajor;
        constexpr static size_t RowAtCompile = Order;
        constexpr static size_t ColAtCompile = Order;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;

        constexpr static bool FastAssign = true;
    };
}

#include "IdentityMatrixImpl/Mul.h"
