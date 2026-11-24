/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "MatrixTrig.h"

namespace Physica {
    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig> || instanceof_tx<M2, MatrixTrig>)
    class GEMM<M1, M2> : public RValueMatrix<GEMM<M1, M2>> {
        using This = GEMM<M1, M2>;
        using Base = RValueMatrix<This>;

        constexpr static bool TrigLHS = instanceof_tx<M1, MatrixTrig>;
        constexpr static bool TrigRHS = instanceof_tx<M2, MatrixTrig>;
        static_assert(!TrigLHS || !TrigRHS, "[Error]: NoImpl");
    protected:
        using typename Base::T;
        using typename Base::Tc;
    private:
        decay_rvalue_t<M1> lhs;
        decay_rvalue_t<M2> rhs;
    public:
        GEMM(M1&& lhs, M2&& rhs);
        GEMM(const This&) = default;
        GEMM(This&&) noexcept = default;
        ~GEMM() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(Matrix auto& target) const;
        void assign_mkl(Matrix auto& target) const noexcept;
        void assign_base(Matrix auto& target) const;
        /* Getters */
        [[nodiscard]] constexpr size_t getRow() const noexcept { return lhs.getRow(); }
        [[nodiscard]] constexpr size_t getCol() const noexcept { return rhs.getCol(); }
        [[nodiscard]] constexpr size_t getOrder() const noexcept { return getRow(); }
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStaticSquare() noexcept;
    };

    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig> || instanceof_tx<M2, MatrixTrig>)
    GEMM<M1, M2>::GEMM(M1&& lhs, M2&& rhs) : lhs(std::forward<M1>(lhs)), rhs(std::forward<M2>(rhs)) {}

    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig> || instanceof_tx<M2, MatrixTrig>)
    void GEMM<M1, M2>::assign(Matrix auto& target) const {
        using Trig = std::conditional_t<TrigLHS, M1, M2>;
        using Data = std::remove_cvref_t<decltype(std::declval<const Trig&>().getExpr())>;
        constexpr bool UseMKL = HasMKL() && !T::isDiffable()
                             && Data::getMinorStrideAtCompile() == 1 && target.getMinorStrideAtCompile() == 1;
        if constexpr (UseMKL)
            assign_mkl(target);
        else
            assign_base(target);
    }

    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig> || instanceof_tx<M2, MatrixTrig>)
    void GEMM<M1, M2>::assign_base(Matrix auto& target) const {
        using Trig = std::conditional_t<TrigLHS, M1, M2>;
        constexpr bool Upper = Traits<Trig>::Upper;
        target.assert_assign(*this);
        target.zeros();

        const size_t maxK = lhs.getCol();
        if constexpr (TrigLHS) {
            for (size_t c = 0; c < getCol(); ++c) {
                auto col = target.col(c);
                for (size_t k = 0; k < maxK; ++k) {
                    if constexpr (Upper)
                        col.head(k + 1) += lhs.col(k).head(k + 1) * rhs.calc(k, c);
                    else
                        col.tail(k) += lhs.col(k).tail(k) * rhs.calc(k, c);
                }
            }
        }
        else {
            for (size_t r = 0; r < getRow(); ++r) {
                auto row = target.row(r);
                for (size_t k = 0; k < maxK; ++k) {
                    if constexpr (Upper)
                        row.tail(k) += lhs.calc(r, k) * rhs.row(k).tail(k);
                    else
                        row.head(k + 1) += lhs.calc(r, k) * rhs.row(k).head(k + 1);
                }
            }
        }
    }

    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig> || instanceof_tx<M2, MatrixTrig>)
    auto&& GEMM<M1, M2>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M1>(self.lhs);
    }

    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig> || instanceof_tx<M2, MatrixTrig>)
    auto&& GEMM<M1, M2>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M2>(self.rhs);
    }

    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig> || instanceof_tx<M2, MatrixTrig>)
    __host__ __device__ consteval bool GEMM<M1, M2>::isStaticSquare() noexcept {
        using LHS = std::remove_cvref_t<M1>;
        using RHS = std::remove_cvref_t<M2>;
        return Base::isStaticSquare() || (LHS::isStaticSquare() && RHS::isStaticSquare());
    }
}

#ifdef PHYSICA_MKL
    #include "GEMM_MKL.h"
#endif
