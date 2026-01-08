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

#include "MatrixTrig.h"

namespace Physica {
    template<Matrix M1, Matrix M2> requires(instanceof<Inverse, M1> && instanceof_tx<MatrixTrig, typename Traits<M1>::ExprType>)
    class GEMM<M1, M2> : public RValueMatrix<GEMM<M1, M2>> {
        using This = GEMM<M1, M2>;
        using Base = RValueMatrix<This>;
    public:
        using Base::isComplex;
    protected:
        using typename Base::T;
    private:
        LazyDestroy<M1&&> inv;
        LazyDestroy<M2&&> rhs;
    public:
        GEMM(M1 inv, M2 rhs);
        GEMM(const This&) = default;
        GEMM(This&&) noexcept = default;
        ~GEMM() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(Matrix auto& target) const;
        void assign_mkl(Matrix auto& target) const noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const { return inv.getRow(); }
        [[nodiscard]] size_t getCol() const { return rhs.getCol(); }
        [[nodiscard]] const auto& getLHS() const noexcept { return inv; }
        [[nodiscard]] const auto& getRHS() const noexcept { return rhs; }
    };

    template<Matrix M1, Matrix M2> requires(instanceof<Inverse, M1> && instanceof_tx<MatrixTrig, typename Traits<M1>::ExprType>)
    GEMM<M1, M2>::GEMM(M1 inv, M2 rhs) : inv(std::forward<M1>(inv)), rhs(std::forward<M2>(rhs)) {}

    template<Matrix M1, Matrix M2> requires(instanceof<Inverse, M1> && instanceof_tx<MatrixTrig, typename Traits<M1>::ExprType>)
    void GEMM<M1, M2>::assign(Matrix auto& target) const {
        if constexpr (HasMKL())
            assign_mkl(target);
        else
            noImpl(__func__);
    }
}

#ifdef PHYSICA_MKL
    #include "InvGEMM_MKL.h"
#endif
