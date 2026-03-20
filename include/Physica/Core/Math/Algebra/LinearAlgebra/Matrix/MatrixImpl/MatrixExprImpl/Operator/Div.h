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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/MatrixExpr.h"

namespace Physica {
    template<class T, class U> requires(Scalar<T> || Scalar<U>)
    class MatrixExpr<ExprID::Div, T, U>
            : public BinaryMatrixExpr<ExprID::Div, T, U> {
        using Base = BinaryMatrixExpr<ExprID::Div, T, U>;
    public:
        using typename Base::ScalarType;
    protected:
        using typename Base::Tv;
    public:
        MatrixExpr(T lhs, U rhs) : Base(std::forward<T>(lhs), std::forward<U>(rhs)) {
            if constexpr (Matrix<T>)
                assert(!Base::getRHS().isZero() && "[Error]: Divide by zero");
        }
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const {
            if constexpr (Matrix<T>)
                return Base::getLHS().calc(row, col) / Base::getRHS();
            else
                return Base::getLHS() / Base::getRHS().calc(row, col);
        }

        [[nodiscard]] Tv calc_value(size_t row, size_t col) const {
            if constexpr (Matrix<T>)
                return Base::getLHS().calc_value(row, col) / Base::getRHS().value();
            else
                return Base::getLHS().value() / Base::getRHS().calc_value(row, col);
        }
    };

    template<class T, class U> requires(Vector<T> || Vector<U>)
    class MatrixExpr<ExprID::Div, T, U>
            : public BinaryMatrixExpr<ExprID::Div, T, U> {
        using Base = BinaryMatrixExpr<ExprID::Div, T, U>;
    public:
        using typename Base::ScalarType;
    protected:
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const {
            if constexpr (Matrix<T>)
                return Base::getLHS().calc(row, col) / Base::getRHS().calc(row);
            else
                return Base::getLHS().calc(row) / Base::getRHS().calc(row, col);
        }

        [[nodiscard]] Tv calc_value(size_t row, size_t col) const {
            if constexpr (Matrix<T>)
                return Base::getLHS().calc_value(row, col) / Base::getRHS().calc_value(row);
            else
                return Base::getLHS().calc_value(row) / Base::getRHS().calc_value(row, col);
        }
    };

    template<Matrix M1, Matrix M2>
    class MatrixExpr<ExprID::Div, M1, M2>
            : public BinaryMatrixExpr<ExprID::Div, M1, M2> {
        using Base = BinaryMatrixExpr<ExprID::Div, M1, M2>;
    public:
        using typename Base::ScalarType;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t row, size_t col) const;
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const;

        void reverse(const Matrix auto& grad) const noexcept;
        using Base::reverse;
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Matrix M1, Matrix M2>
    auto MatrixExpr<ExprID::Div, M1, M2>::calc(size_t row, size_t col) const -> CoDiff<T> {
        assert(!getRHS().calc(row, col).isZero() && "[Error]: Divide by zero");
        return getLHS().calc(row, col) / getRHS().calc(row, col);
    }

    template<Matrix M1, Matrix M2>
    auto MatrixExpr<ExprID::Div, M1, M2>::calc_value(size_t row, size_t col) const -> Tv {
        assert(!getRHS().calc_value(row, col).isZero() && "[Error]: Divide by zero");
        return getLHS().calc_value(row, col) / getRHS().calc_value(row, col);
    }

    template<Matrix M1, Matrix M2>
    void MatrixExpr<ExprID::Div, M1, M2>::reverse(const Matrix auto& grad) const noexcept {
        const auto& lhs = getLHS();
        const auto& rhs = getRHS();
        static_assert(Diffable<M1> && !Diffable<M2>, "[Error]: Not implemented");
        lhs.reverse(divide_elem(grad, rhs));
    }

    template<Matrix M, Scalar U>
    [[nodiscard]] auto divide_elem(M&& m, U&& x) noexcept requires(!DeviceObj<M>) {
        return MatrixExpr<ExprID::Div, M&&, U&&>(std::forward<M>(m), std::forward<U>(x));
    }

    template<Matrix M, Scalar U>
    [[nodiscard]] auto divide_elem(U&& x, M&& m) noexcept requires(!DeviceObj<M>) {
        return MatrixExpr<ExprID::Div, U&&, M&&>(std::forward<U>(x), std::forward<M>(m));
    }

    template<Matrix M, Vector V>
    [[nodiscard]] auto divide_elem(M&& m, V&& x) noexcept requires(!DeviceObj<M> && !DeviceObj<V>) {
        return MatrixExpr<ExprID::Div, M&&, V&&>(std::forward<M>(m), std::forward<V>(x));
    }

    template<Matrix M, Vector V>
    [[nodiscard]] auto divide_elem(V&& x, M&& m) noexcept requires(!DeviceObj<M> && !DeviceObj<V>) {
        return MatrixExpr<ExprID::Div, V&&, M&&>(std::forward<V>(x), std::forward<M>(m));
    }

    template<Matrix M1, Matrix M2>
    [[nodiscard]] auto divide_elem(M1&& m1, M2&& m2) noexcept requires(!DeviceObj<M1> && !DeviceObj<M2>) {
        return MatrixExpr<ExprID::Div, M1&&, M2&&>(std::forward<M1>(m1), std::forward<M2>(m2));
    }
}
