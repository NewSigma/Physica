/*
 * Copyright 2025 Weibo He.
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
    class DiagMatrix : public RValueMatrix<DiagMatrix<T, Order>> {
        using This = DiagMatrix<T, Order>;
        using Base = RValueMatrix<This>;
        using VectorType = DenseVector<T, Order>;
    protected:
        using typename Base::Tv;
    private:
        VectorType diag;
    public:
        DiagMatrix() = default;
        DiagMatrix(VectorType diag_);
        DiagMatrix(const This&) = default;
        DiagMatrix(This&&) noexcept = default;
        ~DiagMatrix() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        template<Vector V>
        [[nodiscard]] V&& operator*(V&& v) const noexcept;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const;
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const;

        [[nodiscard]] const This& transpose() const noexcept { return *this; }
        inline void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getDiag() const noexcept { return diag; }
        [[nodiscard]] size_t getRow() const noexcept;
        [[nodiscard]] size_t getCol() const noexcept;
    };

    template<Scalar T, size_t Order>
    DiagMatrix<T, Order>::DiagMatrix(VectorType diag_) : diag(std::move(diag_)) {}

    template<Scalar T, size_t Order>
    template<Vector V>
    V&& DiagMatrix<T, Order>::operator*(V&& v) const noexcept {
        assert(getCol() == v.getLength());
        return hadamard(diag, v);
    }

    template<Scalar T, size_t Order>
    auto DiagMatrix<T, Order>::calc(size_t row, size_t col) const -> T {
        if (row != col)
            return T(0);
        return diag[row];
    }

    template<Scalar T, size_t Order>
    auto DiagMatrix<T, Order>::calc_value(size_t row, size_t col) const -> Tv {
        return calc(row, col).value();
    }

    template<Scalar T, size_t Order>
    inline void DiagMatrix<T, Order>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(diag, obj.diag);
    }

    template<Scalar T, size_t Order>
    size_t DiagMatrix<T, Order>::getRow() const noexcept {
        if constexpr (Order == Dynamic)
            return diag.getLength();
        return Order;
    }

    template<Scalar T, size_t Order>
    size_t DiagMatrix<T, Order>::getCol() const noexcept {
        return getRow();
    }
}

namespace Physica {
    template<Scalar T, size_t Order>
    class Traits<DiagMatrix<T, Order>> {
    public:
        using ScalarType = T;
        constexpr static int Option = MatrixOption::AnyMajor | MatrixOption::AnyStorage;
        constexpr static size_t RowAtCompile = Order;
        constexpr static size_t ColAtCompile = Order;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;

        constexpr static bool FastAssign = true;
    };
}
