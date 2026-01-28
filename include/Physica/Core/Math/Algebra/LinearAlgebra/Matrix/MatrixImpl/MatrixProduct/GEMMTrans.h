/*
 * Copyright 2025-2026 Weibo He.
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
    template<Matrix M1, Matrix M2> requires(!instanceof<Inverse, M1> && instanceof<Transpose, M2>)
    class GEMM<M1, M2> : public RValueMatrix<GEMM<M1, M2>> {
        using This = GEMM<M1, M2>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
        using typename Base::Tm;
    private:
        LazyDestroy<M1> mat1;
        LazyDestroy<M2> mat2;
    public:
        GEMM(M1&& mat1_, M2&& mat2_);
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
        [[nodiscard]] size_t getCol() const { return mat2.getCol(); }
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
    };

    template<Matrix M1, Matrix M2> requires(!instanceof<Inverse, M1> && instanceof<Transpose, M2>)
    GEMM<M1, M2>::GEMM(M1&& mat1_, M2&& mat2_) : mat1(std::forward<M1>(mat1_)), mat2(std::forward<M2>(mat2_)) {}

    template<Matrix M1, Matrix M2> requires(!instanceof<Inverse, M1> && instanceof<Transpose, M2>)
    void GEMM<M1, M2>::assign(Matrix auto& target) const {
        target.assert_assign(*this);

        using Helper = GEMM<M1, decltype(mat2.getExpr())>;
        if constexpr (Helper::UseMKL(target)) {
            constexpr int Critical = Helper::Critical;
            if (getLHS().getSize() > Critical && getRHS().getSize() > Critical)
                assign_mkl(target);
            else
                Base::assign_base(target);
        }
        else
            Base::assign_base(target);
    }

    template<Matrix M1, Matrix M2> requires(!instanceof<Inverse, M1> && instanceof<Transpose, M2>)
    auto GEMM<M1, M2>::calc(size_t row, size_t col) const -> CoDiff<T> {
        return mat1.row(row) * mat2.col(col);
    }

    template<Matrix M1, Matrix M2> requires(!instanceof<Inverse, M1> && instanceof<Transpose, M2>)
    auto GEMM<M1, M2>::calc_value(size_t row, size_t col) const -> Tv {
        return (getLHS().values() * getRHS().values()).calc(row, col);
    }

    template<Matrix M1, Matrix M2> requires(!instanceof<Inverse, M1> && instanceof<Transpose, M2>)
    auto GEMM<M1, M2>::values() const noexcept {
        return mat1.values() * mat2.values();
    }

    template<Matrix M1, Matrix M2> requires(!instanceof<Inverse, M1> && instanceof<Transpose, M2>)
    auto&& GEMM<M1, M2>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M1>(self.mat1);
    }

    template<Matrix M1, Matrix M2> requires(!instanceof<Inverse, M1> && instanceof<Transpose, M2>)
    auto&& GEMM<M1, M2>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M2>(self.mat2);
    }
}

#ifdef PHYSICA_MKL
    #include "GEMMTrans_MKL.h"
#endif
