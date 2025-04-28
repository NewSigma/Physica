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
    class MatrixProduct<TrigUpper<M1>, M2> : public RValueMatrix<MatrixProduct<TrigUpper<M1>, M2>> {
        using This = MatrixProduct<TrigUpper<M1>, M2>;
        using Base = RValueMatrix<This>;
    public:
        using Base::isComplex;
    protected:
        using typename Base::T;
    private:
        const TrigUpper<M1>& mat1;
        const M2& mat2;
    public:
        MatrixProduct(const TrigUpper<M1>& inv, const M2& mat2_);
        MatrixProduct(const This&) = default;
        MatrixProduct(This&&) noexcept = default;
        ~MatrixProduct() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<Matrix M>
        void assign(M& target) const;
        template<Matrix M>
        void assign_mkl(M& target) const;
        /* Getters */
        [[nodiscard]] size_t getRow() const { return mat1.getRow(); }
        [[nodiscard]] size_t getCol() const { return mat2.getCol(); }
        [[nodiscard]] const auto& getLHS() const noexcept { return mat1; }
        [[nodiscard]] const auto& getRHS() const noexcept { return mat2; }
    };

    template<Matrix M1, Matrix M2>
    MatrixProduct<TrigUpper<M1>, M2>::MatrixProduct(const TrigUpper<M1>& mat1_, const M2& mat2_) : mat1(mat1_), mat2(mat2_) {}

    template<Matrix M1, Matrix M2>
    template<Matrix M>
    void MatrixProduct<TrigUpper<M1>, M2>::assign(M& target) const {
        if constexpr (HasMKL())
            assign_mkl(target);
        else
            noImpl(__func__);
    }
}

#ifdef PHYSICA_MKL
    #include "TrigGEMM_MKL.h"
#endif
