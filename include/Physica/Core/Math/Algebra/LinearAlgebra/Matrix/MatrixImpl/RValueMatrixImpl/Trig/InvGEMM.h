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

#include "MatrixTrig.h"

namespace Physica {
    namespace Internal {
        template<Matrix M>
        consteval bool isInvTrig() noexcept requires(instanceof<M, Inverse> && instanceof_tx<typename Traits<M>::ExprType, MatrixTrig>) {
            return true;
        }

        template<Matrix M>
        consteval bool isInvTrig() noexcept {
            return false;
        }
    }

    template<Matrix M1, Matrix M2> requires(Internal::isInvTrig<M1>() != Internal::isInvTrig<M2>())
    class GEMM<M1, M2> : public RValueMatrix<GEMM<M1, M2>> {
        using This = GEMM<M1, M2>;
        using Base = RValueMatrix<This>;
    public:
        using Base::isComplex;
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
        [[nodiscard]] size_t getRow() const noexcept { return lhs.getRow(); }
        [[nodiscard]] size_t getCol() const noexcept { return rhs.getCol(); }
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStaticSquare() noexcept;
    };

    template<Matrix M1, Matrix M2> requires(Internal::isInvTrig<M1>() != Internal::isInvTrig<M2>())
    GEMM<M1, M2>::GEMM(M1&& lhs, M2&& rhs) : lhs(std::forward<M1>(lhs)), rhs(std::forward<M2>(rhs)) {}

    template<Matrix M1, Matrix M2> requires(Internal::isInvTrig<M1>() != Internal::isInvTrig<M2>())
    void GEMM<M1, M2>::assign(Matrix auto& target) const {
        if constexpr (HasMKL())
            assign_mkl(target);
        else
            assign_base(target);
    }

    template<Matrix M1, Matrix M2> requires(Internal::isInvTrig<M1>() != Internal::isInvTrig<M2>())
    void GEMM<M1, M2>::assign_base(Matrix auto& target) const {
        using M = std::remove_cvref_t<decltype(target)>;
        target.assert_assign(*this);
        if constexpr (Internal::isInvTrig<M1>()) {
            constexpr bool Upper = Traits<M1>::Upper;
            constexpr bool Unit = Traits<M1>::Unit;
            if constexpr (std::same_as<std::remove_cvref_t<M2>, M>) {
                if (&rhs != &target)
                    rhs.assign(target);
            }
            else
                rhs.assign(target);

            const auto& trig = lhs.getExpr();
            const size_t order = getRow();
            if constexpr (Upper) {
                for (size_t i = order; i > 0; --i) {
                    const size_t j = i - 1;
                    auto row = target.row(j);
                    if constexpr (!Unit)
                        row /= trig.calc(j, j);
                    if (j > 0)
                        target.topRows(j) -= trig.col(j).head(j) * row.transpose();
                }
            }
            else {
                for (size_t i = 0; i < order; ++i) {
                    auto row = target.row(i);
                    if constexpr (!Unit)
                        row /= trig.calc(i, i);
                    if (i + 1 < order)
                        target.bottomRows(i + 1) -= trig.col(i).tail(i + 1) * row.transpose();
                }
            }
        }
        else {
            constexpr bool Upper = Traits<M2>::Upper;
            constexpr bool Unit = Traits<M2>::Unit;
            if constexpr (std::same_as<std::remove_cvref_t<M1>, M>) {
                if (&lhs != &target)
                    lhs.assign(target);
            }
            else
                lhs.assign(target);

            const auto& trig = rhs.getExpr();
            const size_t order = getCol();
            if constexpr (Upper) {
                for (size_t i = 0; i < order; ++i) {
                    auto col = target.col(i);
                    if constexpr (!Unit)
                        col /= trig.calc(i, i);
                    if (i + 1 < order)
                        target.rightCols(i + 1) -= col * trig.row(i).tail(i + 1).transpose();
                }
            }
            else {
                for (size_t i = order; i > 0; --i) {
                    const size_t j = i - 1;
                    auto col = target.col(j);
                    if constexpr (!Unit)
                        col /= trig.calc(j, j);
                    if (j > 0)
                        target.leftCols(j) -= col * trig.row(j).head(j).transpose();
                }
            }
        }
    }

    template<Matrix M1, Matrix M2> requires(Internal::isInvTrig<M1>() != Internal::isInvTrig<M2>())
    auto&& GEMM<M1, M2>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M1>(self.lhs);
    }

    template<Matrix M1, Matrix M2> requires(Internal::isInvTrig<M1>() != Internal::isInvTrig<M2>())
    auto&& GEMM<M1, M2>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M2>(self.rhs);
    }

    template<Matrix M1, Matrix M2> requires(Internal::isInvTrig<M1>() != Internal::isInvTrig<M2>())
    __host__ __device__ consteval bool GEMM<M1, M2>::isStaticSquare() noexcept {
        using LHS = std::remove_cvref_t<M1>;
        using RHS = std::remove_cvref_t<M2>;
        return Base::isStaticSquare() || (LHS::isStaticSquare() && RHS::isStaticSquare());
    }
}

#ifdef PHYSICA_MKL
    #include "InvGEMM_MKL.h"
#endif
