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

#include "../DiagMatrix.h"

namespace Physica {
    template<Matrix M1, Scalar U, size_t Order>
    class GEMM<M1, DiagMatrix<U, Order>> : public RValueMatrix<GEMM<M1, DiagMatrix<U, Order>>> {
        using This = GEMM<M1, DiagMatrix<U, Order>>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
    private:
        const M1& lhs;
        const DiagMatrix<U, Order>& rhs;
    public:
        GEMM(const M1& lhs_, const DiagMatrix<U, Order>& rhs_);
        GEMM(const This&) = default;
        GEMM(This&&) noexcept = default;
        ~GEMM() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<Matrix M>
        void assign(M& target) const;

        [[nodiscard]] T calc(size_t row, size_t col) const;
        [[nodiscard]] T calc_value(size_t row, size_t col) const;
        /* Getters */
        [[nodiscard]] size_t getRow() const { return lhs.getRow(); }
        [[nodiscard]] size_t getCol() const { return rhs.getCol(); }
        [[nodiscard]] const M1& getLHS() const noexcept { return lhs; }
        [[nodiscard]] const auto& getRHS() const noexcept { return rhs; }
    };

    template<Matrix M1, Scalar U, size_t Order>
    GEMM<M1, DiagMatrix<U, Order>>::GEMM(
            const M1& lhs_, const DiagMatrix<U, Order>& rhs_) : lhs(lhs_), rhs(rhs_) {}

    template<Matrix M1, Scalar U, size_t Order>
    template<Matrix M>
    void GEMM<M1, DiagMatrix<U, Order>>::assign(M& target) const {
        for (size_t i = 0; i < getCol(); ++i)
            target.col(i) = lhs.col(i) * rhs.diag()[i];
    }

    template<Matrix M1, Scalar U, size_t Order>
    auto GEMM<M1, DiagMatrix<U, Order>>::calc(size_t row, size_t col) const -> T {
        return lhs.calc(row, col) * rhs.diag()[col];
    }

    template<Matrix M1, Scalar U, size_t Order>
    auto GEMM<M1, DiagMatrix<U, Order>>::calc_value(size_t row, size_t col) const -> T {
        return lhs.calc_value(row, col) * rhs.diag()[col].value();
    }

    template<Scalar U, size_t Order, Matrix M2>
    class GEMM<DiagMatrix<U, Order>, M2> : public RValueMatrix<GEMM<DiagMatrix<U, Order>, M2>> {
        using This = GEMM<DiagMatrix<U, Order>, M2>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
    private:
        const DiagMatrix<U, Order>& lhs;
        const M2& rhs;
    public:
        GEMM(const DiagMatrix<U, Order>& lhs_, const M2& rhs_);
        GEMM(const This&) = default;
        GEMM(This&&) noexcept = default;
        ~GEMM() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<Matrix M>
        void assign(M& target) const;

        [[nodiscard]] T calc(size_t row, size_t col) const;
        [[nodiscard]] T calc_value(size_t row, size_t col) const;
        /* Getters */
        [[nodiscard]] size_t getRow() const { return lhs.getRow(); }
        [[nodiscard]] size_t getCol() const { return rhs.getCol(); }
        [[nodiscard]] const auto& getLHS() const noexcept { return lhs; }
        [[nodiscard]] const auto& getRHS() const noexcept { return rhs; }
    };

    template<Scalar U, size_t Order, Matrix M2>
    GEMM<DiagMatrix<U, Order>, M2>::GEMM(
            const DiagMatrix<U, Order>& lhs_, const M2& rhs_) : lhs(lhs_), rhs(rhs_) {}

    template<Scalar U, size_t Order, Matrix M2>
    template<Matrix M>
    void GEMM<DiagMatrix<U, Order>, M2>::assign(M& target) const {
        for (size_t i = 0; i < getRow(); ++i)
            target.row(i) = lhs.diag()[i] * rhs.row(i);
    }

    template<Scalar U, size_t Order, Matrix M2>
    auto GEMM<DiagMatrix<U, Order>, M2>::calc(size_t row, size_t col) const -> T {
        return lhs.diag()[row] * rhs.calc(row, col);
    }

    template<Scalar U, size_t Order, Matrix M2>
    auto GEMM<DiagMatrix<U, Order>, M2>::calc_value(size_t row, size_t col) const -> T {
        return lhs.diag()[row].value() * rhs.calc_value(row, col);
    }
}
