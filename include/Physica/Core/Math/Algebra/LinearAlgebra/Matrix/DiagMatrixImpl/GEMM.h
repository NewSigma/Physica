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
    template<Matrix M1, Matrix M2> requires(instanceof_tx<DiagMatrix, M2>)
    class GEMM<M1, M2> : public RValueMatrix<GEMM<M1, M2>> {
        using This = GEMM<M1, M2>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
    private:
        LazyDestroy<M1> lhs;
        LazyDestroy<M2> rhs;
    public:
        GEMM(M1&& lhs_, M2&& rhs_);
        GEMM(const This&) = default;
        GEMM(This&&) noexcept = default;
        ~GEMM() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(Matrix auto& target) const;

        [[nodiscard]] T calc(size_t row, size_t col) const;

        [[nodiscard]] auto values(this auto&& self) noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const { return lhs.getRow(); }
        [[nodiscard]] size_t getCol() const { return rhs.getCol(); }
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept;
    };

    template<Matrix M1, Matrix M2> requires(instanceof_tx<DiagMatrix, M2>)
    GEMM<M1, M2>::GEMM(M1&& lhs_, M2&& rhs_) : lhs(std::forward<M1>(lhs_)), rhs(std::forward<M2>(rhs_)) {}

    template<Matrix M1, Matrix M2> requires(instanceof_tx<DiagMatrix, M2>)
    void GEMM<M1, M2>::assign(Matrix auto& target) const {
        for (size_t i = 0; i < getCol(); ++i)
            target.col(i) = lhs.col(i) * rhs.diag()[i];
    }

    template<Matrix M1, Matrix M2> requires(instanceof_tx<DiagMatrix, M2>)
    auto GEMM<M1, M2>::calc(size_t row, size_t col) const -> T {
        return lhs.calc(row, col) * rhs.diag()[col];
    }

    template<Matrix M1, Matrix M2> requires(instanceof_tx<DiagMatrix, M2>)
    auto GEMM<M1, M2>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS() * std::forward<Self>(self).getRHS();
    }

    template<Matrix M1, Matrix M2> requires(instanceof_tx<DiagMatrix, M2>)
    auto&& GEMM<M1, M2>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M1>(self.lhs);
    }

    template<Matrix M1, Matrix M2> requires(instanceof_tx<DiagMatrix, M2>)
    auto&& GEMM<M1, M2>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M2>(self.rhs);
    }

    template<Matrix M1, Matrix M2> requires(instanceof_tx<DiagMatrix, M2>)
    __host__ __device__ consteval int GEMM<M1, M2>::getMajor() noexcept {
        return MatrixMajor::getMajor<M1>();
    }

    template<Matrix M1, Matrix M2> requires(instanceof_tx<DiagMatrix, M1>)
    class GEMM<M1, M2> : public RValueMatrix<GEMM<M1, M2>> {
        using This = GEMM<M1, M2>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
    private:
        LazyDestroy<M1> lhs;
        LazyDestroy<M2> rhs;
    public:
        GEMM(M1&& lhs_, M2&& rhs_);
        GEMM(const This&) = default;
        GEMM(This&&) noexcept = default;
        ~GEMM() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(Matrix auto& target) const;

        [[nodiscard]] T calc(size_t row, size_t col) const;

        [[nodiscard]] auto values(this auto&& self) noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const { return lhs.getRow(); }
        [[nodiscard]] size_t getCol() const { return rhs.getCol(); }
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept;
    };

    template<Matrix M1, Matrix M2> requires(instanceof_tx<DiagMatrix, M1>)
    GEMM<M1, M2>::GEMM(M1&& lhs_, M2&& rhs_) : lhs(std::forward<M1>(lhs_)), rhs(std::forward<M2>(rhs_)) {}

    template<Matrix M1, Matrix M2> requires(instanceof_tx<DiagMatrix, M1>)
    void GEMM<M1, M2>::assign(Matrix auto& target) const {
        for (size_t i = 0; i < getRow(); ++i)
            target.row(i) = lhs.diag()[i] * rhs.row(i);
    }

    template<Matrix M1, Matrix M2> requires(instanceof_tx<DiagMatrix, M1>)
    auto GEMM<M1, M2>::calc(size_t row, size_t col) const -> T {
        return lhs.diag()[row] * rhs.calc(row, col);
    }

    template<Matrix M1, Matrix M2> requires(instanceof_tx<DiagMatrix, M1>)
    auto GEMM<M1, M2>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS() * std::forward<Self>(self).getRHS();
    }

    template<Matrix M1, Matrix M2> requires(instanceof_tx<DiagMatrix, M1>)
    auto&& GEMM<M1, M2>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M1>(self.lhs);
    }

    template<Matrix M1, Matrix M2> requires(instanceof_tx<DiagMatrix, M1>)
    auto&& GEMM<M1, M2>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M2>(self.rhs);
    }

    template<Matrix M1, Matrix M2> requires(instanceof_tx<DiagMatrix, M1>)
    __host__ __device__ consteval int GEMM<M1, M2>::getMajor() noexcept {
        return MatrixMajor::getMajor<M2>();
    }
}
