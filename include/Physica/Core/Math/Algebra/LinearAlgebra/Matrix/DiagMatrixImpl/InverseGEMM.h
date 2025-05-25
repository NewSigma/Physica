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
    class GEMM<M1, Inverse<DiagMatrix<U, Order>>> : public RValueMatrix<GEMM<M1, Inverse<DiagMatrix<U, Order>>>> {
        using This = GEMM<M1, Inverse<DiagMatrix<U, Order>>>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
    private:
        const M1& mat1;
        const DiagMatrix<U, Order>& mat2;
    public:
        GEMM(const M1& inv, const Inverse<DiagMatrix<U, Order>>& mat2_);
        GEMM(const This&) = default;
        GEMM(This&&) noexcept = default;
        ~GEMM() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<Matrix M>
        void assign(M& target) const;
        /* Getters */
        [[nodiscard]] size_t getRow() const { return mat1.getRow(); }
        [[nodiscard]] size_t getCol() const { return mat2.getCol(); }
        [[nodiscard]] const M1& getLHS() const noexcept { return mat1; }
        [[nodiscard]] auto getRHS() const noexcept { return mat2.inverse(); }
    };

    template<Matrix M1, Scalar U, size_t Order>
    GEMM<M1, Inverse<DiagMatrix<U, Order>>>::GEMM(
            const M1& mat1_, const Inverse<DiagMatrix<U, Order>>& inv) : mat1(mat1_), mat2(inv.getExpr()) {}

    template<Matrix M1, Scalar U, size_t Order>
    template<Matrix M>
    void GEMM<M1, Inverse<DiagMatrix<U, Order>>>::assign(M& target) const {
        for (size_t i = 0; i < getCol(); ++i)
            target.col(i) = mat1.col(i) * reciprocal(mat2.diag()[i]);
    }
}
