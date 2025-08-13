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

#include "Trig.h"

namespace Physica {
    template<Matrix M1, Matrix M2>
    class GEMM<TrigUpper<M1>, M2> : public RValueMatrix<GEMM<TrigUpper<M1>, M2>> {
        using This = GEMM<TrigUpper<M1>, M2>;
        using Base = RValueMatrix<This>;
    public:
        using Base::isComplex;
    protected:
        using typename Base::T;
    private:
        const TrigUpper<M1>& mat1;
        const M2& mat2;
    public:
        GEMM(const TrigUpper<M1>& inv, const M2& mat2_);
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
        [[nodiscard]] size_t getRow() const { return mat1.getRow(); }
        [[nodiscard]] size_t getCol() const { return mat2.getCol(); }
        [[nodiscard]] const auto& getLHS() const noexcept { return mat1; }
        [[nodiscard]] const auto& getRHS() const noexcept { return mat2; }
    };

    template<Matrix M1, Matrix M2>
    GEMM<TrigUpper<M1>, M2>::GEMM(const TrigUpper<M1>& mat1_, const M2& mat2_) : mat1(mat1_), mat2(mat2_) {}

    template<Matrix M1, Matrix M2>
    void GEMM<TrigUpper<M1>, M2>::assign(Matrix auto& target) const {
        if constexpr (HasMKL())
            assign_mkl(target);
        else
            noImpl(__func__);
    }
}

#ifdef PHYSICA_MKL
    #include "TrigGEMM_MKL.h"
#endif
