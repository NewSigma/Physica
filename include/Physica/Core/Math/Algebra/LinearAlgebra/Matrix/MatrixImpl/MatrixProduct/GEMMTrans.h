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

#include "GEMM.h"

namespace Physica {
    template<Matrix M1, Matrix M2> requires(!instanceof<Inverse, M1>)
    class GEMM<M1, Transpose<M2>> : public RValueMatrix<GEMM<M1, Transpose<M2>>> {
        using This = GEMM<M1, Transpose<M2>>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
        using typename Base::Tm;
    private:
        const M1& mat1;
        const M2& mat2;
    public:
        GEMM(const M1& mat1_, const Transpose<M2>& mat2_);
        GEMM(const This&) = default;
        GEMM(This&&) noexcept = default;
        ~GEMM() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(Matrix auto& target) const;
        void assign_mkl(Matrix auto& target) const noexcept;

        [[nodiscard]] CoDiff<T> calc(size_t row, size_t col) const;
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const;

        [[nodiscard]] auto values() const noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const { return mat1.getRow(); }
        [[nodiscard]] size_t getCol() const { return mat2.getRow(); }
        [[nodiscard]] const auto& getLHS() const noexcept { return mat1; }
        [[nodiscard]] decltype(auto) getRHS() const noexcept { return mat2.transpose(); }
    };

    template<Matrix M1, Matrix M2> requires(!instanceof<Inverse, M1>)
    GEMM<M1, Transpose<M2>>::GEMM(const M1& mat1_, const Transpose<M2>& mat2_) : mat1(mat1_), mat2(mat2_.getExpr()) {}

    template<Matrix M1, Matrix M2> requires(!instanceof<Inverse, M1>)
    void GEMM<M1, Transpose<M2>>::assign(Matrix auto& target) const {
        target.assert_assign(*this);
        if constexpr (GEMM<M1, M2>::template UseMKL<decltype(target)>()) {
            constexpr int Critical = GEMM<M1, M2>::Critical;
            if (getLHS().getSize() > Critical && getRHS().getSize() > Critical)
                assign_mkl(target);
            else
                Base::assign_base(target);
        }
        else
            Base::assign_base(target);
    }

    template<Matrix M1, Matrix M2> requires(!instanceof<Inverse, M1>)
    auto GEMM<M1, Transpose<M2>>::calc(size_t row, size_t col) const -> CoDiff<T> {
        return mat1.row(row) * mat2.row(col);
    }

    template<Matrix M1, Matrix M2> requires(!instanceof<Inverse, M1>)
    auto GEMM<M1, Transpose<M2>>::calc_value(size_t row, size_t col) const -> Tv {
        return (getLHS().values() * getRHS().values()).calc(row, col);
    }

    template<Matrix M1, Matrix M2> requires(!instanceof<Inverse, M1>)
    auto GEMM<M1, Transpose<M2>>::values() const noexcept {
        return mat1 * mat2.transpose();
    }
}

#ifdef PHYSICA_MKL
    #include "GEMMTrans_MKL.h"
#endif
