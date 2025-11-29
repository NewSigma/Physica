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

#include "../MatrixExpr.h"

namespace Physica {
    template<class U, class V>
    class MatrixExpr<ExprID::Sub, U, V>
            : public BinaryMatrixExpr<ExprID::Sub, U, V> {
        static_assert(Scalar<U> || Scalar<V>, "[Error]: Either types should be Scalar");

        using Base = BinaryMatrixExpr<ExprID::Sub, U, V>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        void assign(Matrix auto& target) const;

        [[nodiscard]] T calc(size_t row, size_t col) const;
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const;
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<class U, class V>
    void MatrixExpr<ExprID::Sub, U, V>::assign(Matrix auto& target) const {
        if constexpr (Matrix<U>) {
            if constexpr (MatrixOption::isSameMajor<U, decltype(target)>())
                (getLHS().flatten() - getRHS()).assign(target.flatten());
            else
                Base::assign(target);
        }
        else {
            if constexpr (MatrixOption::isSameMajor<V, decltype(target)>())
                (getLHS() - getRHS().flatten()).assign(target.flatten());
            else
                Base::assign(target);
        }
    }

    template<class U, class V>
    auto MatrixExpr<ExprID::Sub, U, V>::calc(size_t row, size_t col) const -> T {
        if constexpr (Matrix<U>)
            return Base::getLHS().calc(row, col) - Base::getRHS();
        else
            return Base::getLHS() - Base::getRHS().calc(row, col);
    }

    template<class U, class V>
    auto MatrixExpr<ExprID::Sub, U, V>::calc_value(size_t row, size_t col) const -> Tv {
        if constexpr (Matrix<U>)
            return Base::getLHS().calc_value(row, col) - Base::getRHS().value();
        else
            return Base::getLHS().value() - Base::getRHS().calc_value(row, col);
    }

    template<Matrix M1, Matrix M2>
    class MatrixExpr<ExprID::Sub, M1, M2>
            : public BinaryMatrixExpr<ExprID::Sub, M1, M2> {
        using Base = BinaryMatrixExpr<ExprID::Sub, M1, M2>;
        using This = MatrixExpr<ExprID::Sub, M1, M2>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        void assign(Matrix auto& target) const;

        [[nodiscard]] CoDiff<T> calc(size_t row, size_t col) const;
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const;

        using Base::reverse;
        void reverse(const auto& grad) const noexcept;

        [[nodiscard]] auto values() const noexcept;
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Matrix M1, Matrix M2>
    void MatrixExpr<ExprID::Sub, M1, M2>::assign(Matrix auto& target) const {
        constexpr bool SameMajor1 = MatrixOption::isSameMajor<M1, decltype(target)>();
        constexpr bool SameMajor2 = MatrixOption::isSameMajor<M2, decltype(target)>();
        if constexpr (SameMajor1 && SameMajor2)
            (getLHS().flatten() - getRHS().flatten()).assign(target.flatten());
        else
            Base::assign(target);
    }

    template<Matrix M1, Matrix M2>
    auto MatrixExpr<ExprID::Sub, M1, M2>::calc(size_t row, size_t col) const -> CoDiff<T> {
        return getLHS().calc(row, col) - getRHS().calc(row, col);
    }

    template<Matrix M1, Matrix M2>
    auto MatrixExpr<ExprID::Sub, M1, M2>::calc_value(size_t row, size_t col) const -> Tv {
        return getLHS().calc_value(row, col) - getRHS().calc_value(row, col);
    }

    template<Matrix M1, Matrix M2>
    void MatrixExpr<ExprID::Sub, M1, M2>::reverse(const auto& grad) const noexcept {
        if constexpr (ReverseDiff<M1>)
            getLHS().reverse(grad);
        if constexpr (ReverseDiff<M2>)
            getRHS().reverse(-grad);
    }

    template<Matrix M1, Matrix M2>
    [[nodiscard]] auto MatrixExpr<ExprID::Sub, M1, M2>::values() const noexcept {
        return getLHS().values() - getRHS().values();
    }

    template<Matrix M, Scalar U>
    [[nodiscard]] auto operator-(M&& m, U&& x) noexcept requires(!CUDA<M>) {
        return MatrixExpr<ExprID::Sub, M&&, U&&>(std::forward<M>(m), std::forward<U>(x));
    }

    template<Matrix M, Scalar U>
    [[nodiscard]] auto operator-(U&& x, M&& m) noexcept requires(!CUDA<M>) {
        return MatrixExpr<ExprID::Sub, U&&, M&&>(std::forward<U>(x), std::forward<M>(m));
    }

    template<Matrix M1, Matrix M2>
    [[nodiscard]] auto operator-(M1&& m1, M2&& m2) noexcept requires(!CUDA<M1> && !CUDA<M2>) {
        return MatrixExpr<ExprID::Sub, M1&&, M2&&>(std::forward<M1>(m1), std::forward<M2>(m2));
    }
}
