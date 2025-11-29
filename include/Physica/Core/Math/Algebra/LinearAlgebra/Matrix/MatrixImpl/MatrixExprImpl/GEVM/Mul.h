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

#include "../Mul.h"

namespace Physica {
    template<Matrix M, Scalar U> requires(instanceof<GEVM, M>)
    class MatrixExpr<ExprID::Mul, M, U> : public BinaryMatrixExpr<ExprID::Mul, M, U> {
        using Base = BinaryMatrixExpr<ExprID::Mul, M, U>;
        using This = MatrixExpr<ExprID::Mul, M, U>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        void assign(Matrix auto& target) const;
        void assign_mkl(Matrix auto& target) const;
        void assign_base(Matrix auto& target) const;
        void assign_add(Matrix auto& target) const;

        [[nodiscard]] CoDiff<T> calc(size_t row, size_t col) const;
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const;
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Matrix M, Scalar U> requires(instanceof<GEVM, M>)
    void MatrixExpr<ExprID::Mul, M, U>::assign(Matrix auto& target) const {
        assign_base(target);
    }

    template<Matrix M, Scalar U> requires(instanceof<GEVM, M>)
    void MatrixExpr<ExprID::Mul, M, U>::assign_base(Matrix auto& target) const {
        const auto& vec = getLHS().getLHS();
        const auto& mat = getLHS().getRHS();
        ((getRHS() * vec) * mat).assign(target);
    }

    template<Matrix M, Scalar U> requires(instanceof<GEVM, M>)
    void MatrixExpr<ExprID::Mul, M, U>::assign_add(Matrix auto& target) const {
        const auto& vec = getLHS().getLHS();
        const auto& mat = getLHS().getRHS();
        ((getRHS() * vec) * mat).assign_add(target);
    }

    template<Matrix M, Scalar U> requires(instanceof<GEVM, M>)
    auto MatrixExpr<ExprID::Mul, M, U>::calc(size_t row, size_t col) const -> CoDiff<T> {
        return getLHS().calc(row, col) * getRHS();
    }

    template<Matrix M, Scalar U> requires(instanceof<GEVM, M>)
    auto MatrixExpr<ExprID::Mul, M, U>::calc_value(size_t row, size_t col) const -> Tv {
        return getLHS().value(row, col) * getRHS().value();
    }
}
