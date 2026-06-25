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

#include "../DiagMatrix.h"

namespace Physica {
    template<Matrix M1, Matrix M2> requires(instanceof<M2, Inverse> && instanceof_tx<typename Traits<M2>::ExprType, DiagMatrix>)
    class GEMM<M1, M2> : public RValueMatrix<GEMM<M1, M2>> {
        using This = GEMM<M1, M2>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
    private:
        decay_rvalue_t<M1> mat1;
        decay_rvalue_t<M2> mat2;
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
        /* Getters */
        [[nodiscard]] size_t getRow() const { return mat1.getRow(); }
        [[nodiscard]] size_t getCol() const { return mat2.getCol(); }
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStaticSquare() noexcept;
    };

    template<Matrix M1, Matrix M2> requires(instanceof<M2, Inverse> && instanceof_tx<typename Traits<M2>::ExprType, DiagMatrix>)
    GEMM<M1, M2>::GEMM(M1&& mat1_, M2&& mat2_) : mat1(std::forward<M1>(mat1_)), mat2(std::forward<M2>(mat2_)) {}

    template<Matrix M1, Matrix M2> requires(instanceof<M2, Inverse> && instanceof_tx<typename Traits<M2>::ExprType, DiagMatrix>)
    void GEMM<M1, M2>::assign(Matrix auto& target) const {
        for (size_t i = 0; i < getCol(); ++i)
            target.col(i) = mat1.col(i) * reciprocal(mat2.diag()[i]);
    }

    template<Matrix M1, Matrix M2> requires(instanceof<M2, Inverse> && instanceof_tx<typename Traits<M2>::ExprType, DiagMatrix>)
    auto&& GEMM<M1, M2>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M1>(self.mat1);
    }

    template<Matrix M1, Matrix M2> requires(instanceof<M2, Inverse> && instanceof_tx<typename Traits<M2>::ExprType, DiagMatrix>)
    auto&& GEMM<M1, M2>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M2>(self.mat2);
    }

    template<Matrix M1, Matrix M2> requires(instanceof<M2, Inverse> && instanceof_tx<typename Traits<M2>::ExprType, DiagMatrix>)
    __host__ __device__ consteval bool GEMM<M1, M2>::isStaticSquare() noexcept {
        using LHS = std::remove_cvref_t<M1>;
        using RHS = std::remove_cvref_t<M2>;
        return Base::isStaticSquare() || (LHS::isStaticSquare() && RHS::isStaticSquare());
    }

    template<Matrix M1, Matrix M2> requires(instanceof<M1, Inverse> && instanceof_tx<typename Traits<M1>::ExprType, DiagMatrix>)
    class GEMM<M1, M2> : public RValueMatrix<GEMM<M1, M2>> {
        using This = GEMM<M1, M2>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
    private:
        decay_rvalue_t<M1> mat1;
        decay_rvalue_t<M2> mat2;
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
        /* Getters */
        [[nodiscard]] size_t getRow() const { return mat1.getRow(); }
        [[nodiscard]] size_t getCol() const { return mat2.getCol(); }
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStaticSquare() noexcept;
    };

    template<Matrix M1, Matrix M2> requires(instanceof<M1, Inverse> && instanceof_tx<typename Traits<M1>::ExprType, DiagMatrix>)
    GEMM<M1, M2>::GEMM(M1&& mat1_, M2&& mat2_) : mat1(std::forward<M1>(mat1_)), mat2(std::forward<M2>(mat2_)) {}

    template<Matrix M1, Matrix M2> requires(instanceof<M1, Inverse> && instanceof_tx<typename Traits<M1>::ExprType, DiagMatrix>)
    void GEMM<M1, M2>::assign(Matrix auto& target) const {
        for (size_t i = 0; i < getCol(); ++i)
            target.row(i) = reciprocal(mat1.diag()[i]) * mat2.row(i);
    }

    template<Matrix M1, Matrix M2> requires(instanceof<M1, Inverse> && instanceof_tx<typename Traits<M1>::ExprType, DiagMatrix>)
    auto&& GEMM<M1, M2>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M1>(self.mat1);
    }

    template<Matrix M1, Matrix M2> requires(instanceof<M1, Inverse> && instanceof_tx<typename Traits<M1>::ExprType, DiagMatrix>)
    auto&& GEMM<M1, M2>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M2>(self.mat2);
    }

    template<Matrix M1, Matrix M2> requires(instanceof<M1, Inverse> && instanceof_tx<typename Traits<M1>::ExprType, DiagMatrix>)
    __host__ __device__ consteval bool GEMM<M1, M2>::isStaticSquare() noexcept {
        using LHS = std::remove_cvref_t<M1>;
        using RHS = std::remove_cvref_t<M2>;
        return Base::isStaticSquare() || (LHS::isStaticSquare() && RHS::isStaticSquare());
    }
}
