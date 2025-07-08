/*
 * Copyright 2024-2025 Weibo He.
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
    class UnitMatrix : public RValueMatrix<UnitMatrix<T, Order>> {
        using This = UnitMatrix<T, Order>;
        using Base = RValueMatrix<This>;
        using IndexType = std::conditional<Order == Dynamic, size_t, PlainStruct<void>>::type;
    protected:
        using typename Base::Tv;
    private:
        [[no_unique_address]] IndexType order;
    public:
        UnitMatrix() = default;
        UnitMatrix(size_t order_);
        UnitMatrix(const This&) = default;
        UnitMatrix(This&&) noexcept = default;
        ~UnitMatrix() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        template<Vector V>
        [[nodiscard]] V&& operator*(V&& v) const noexcept;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const { return T(row == col ? 1 : 0); }
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const { return Tv(row == col ? 1 : 0); }

        [[nodiscard]] const This& transpose() const noexcept { return *this; }
        [[nodiscard]] const This& conjugate() const noexcept { return *this; }
        [[nodiscard]] const This& hermite() const noexcept { return *this; }
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept;
        [[nodiscard]] size_t getCol() const noexcept;
    };

    template<Scalar T, size_t Order>
    UnitMatrix<T, Order>::UnitMatrix(size_t order_) : order(order_) {
        assert((Order == Dynamic || Order == order_) && "[Error]: tparam and param is not consistent");
    }

    template<Scalar T, size_t Order>
    template<Vector V>
    V&& UnitMatrix<T, Order>::operator*(V&& v) const noexcept {
        assert(getCol() == v.getLength());
        return std::forward<V>(v);
    }

    template<Scalar T, size_t Order>
    void UnitMatrix<T, Order>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(order, obj.order);
    }

    template<Scalar T, size_t Order>
    size_t UnitMatrix<T, Order>::getRow() const noexcept {
        if constexpr (Order == Dynamic)
            return order;
        return Order;
    }

    template<Scalar T, size_t Order>
    size_t UnitMatrix<T, Order>::getCol() const noexcept {
        return getRow();
    }
}

namespace Physica {
    template<Scalar T, size_t Order>
    class Traits<UnitMatrix<T, Order>> {
    public:
        using ScalarType = T;
        constexpr static int Option = MatrixOption::AnyMajor | MatrixOption::AnyStorage;
        constexpr static size_t RowAtCompile = Order;
        constexpr static size_t ColAtCompile = Order;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;

        constexpr static bool FastAssign = true;
    };
}
