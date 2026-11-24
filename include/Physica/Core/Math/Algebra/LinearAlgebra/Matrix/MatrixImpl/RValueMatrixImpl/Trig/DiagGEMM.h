/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "MatrixTrig.h"

namespace Physica {
    template<Scalar T, size_t Order> class DiagMatrix;

    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig> && instanceof_tx<M2, DiagMatrix>)
    class GEMM<M1, M2> : public RValueMatrix<GEMM<M1, M2>> {
        using This = GEMM<M1, M2>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
    private:
        decay_rvalue_t<M1> trig;
        decay_rvalue_t<M2> diag;
    public:
        GEMM(M1&& trig, M2&& diag);
        GEMM(const This&) = default;
        GEMM(This&&) noexcept = default;
        ~GEMM() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(Matrix auto& target) const;
        /* Getters */
        [[nodiscard]] constexpr size_t getRow() const noexcept { return trig.getRow(); }
        [[nodiscard]] constexpr size_t getCol() const noexcept { return diag.getCol(); }
        [[nodiscard]] constexpr size_t getOrder() const noexcept { return diag.getOrder(); }
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStaticSquare() noexcept;
    };

    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig> && instanceof_tx<M2, DiagMatrix>)
    GEMM<M1, M2>::GEMM(M1&& trig, M2&& diag) : trig(std::forward<M1>(trig)), diag(std::forward<M2>(diag)) {}

    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig> && instanceof_tx<M2, DiagMatrix>)
    void GEMM<M1, M2>::assign(Matrix auto& target) const {
        constexpr bool Upper = Traits<M1>::Upper;
        target.zeros();
        for (size_t i = 0; i < getCol(); ++i) {
            if constexpr (Upper)
                target.col(i).head(i + 1) = trig.col(i).head(i + 1) * diag.diag()[i];
            else
                target.col(i).tail(i) = trig.col(i).tail(i) * diag.diag()[i];
        }
    }

    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig> && instanceof_tx<M2, DiagMatrix>)
    auto&& GEMM<M1, M2>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M1>(self.trig);
    }

    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig> && instanceof_tx<M2, DiagMatrix>)
    auto&& GEMM<M1, M2>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M2>(self.diag);
    }

    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig> && instanceof_tx<M2, DiagMatrix>)
    __host__ __device__ consteval bool GEMM<M1, M2>::isStaticSquare() noexcept {
        return std::remove_cvref_t<M1>::isStaticSquare();
    }
}
