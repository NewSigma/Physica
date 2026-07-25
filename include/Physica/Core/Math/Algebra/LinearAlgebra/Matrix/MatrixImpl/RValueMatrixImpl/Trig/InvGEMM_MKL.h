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

#include "InvGEMM.h"

namespace Physica {
    template<Matrix M1, Matrix M2> requires(Internal::isInvTrig<M1>() != Internal::isInvTrig<M2>())
    void GEMM<M1, M2>::assign_mkl(Matrix auto& target) const noexcept {
        using M = std::remove_cvref_t<decltype(target)>;
        using Tm = decltype(std::declval<T>().toMKL());
        if constexpr (Internal::isInvTrig<M1>()) {
            if constexpr (std::same_as<std::remove_cvref_t<M2>, M>) {
                if (&rhs != &target)
                    rhs.assign(target);
            }
            else
                rhs.assign(target);
        }
        else {
            if constexpr (std::same_as<std::remove_cvref_t<M1>, M>) {
                if (&lhs != &target)
                    lhs.assign(target);
            }
            else
                lhs.assign(target);
        }

        constexpr auto Layout = MatrixMajor::isRowMatrix<M>() ? CblasRowMajor : CblasColMajor;
        constexpr auto Side = Internal::isInvTrig<M1>() ? CblasLeft : CblasRight;
        constexpr auto TransA = []() consteval noexcept {
            if constexpr (Internal::isInvTrig<M1>())
                return MatrixMajor::isSameMajor<M1, M>() ? CblasNoTrans : CblasTrans;
            else
                return MatrixMajor::isSameMajor<M2, M>() ? CblasNoTrans : CblasTrans;
        }();

        constexpr auto UploNoTrans = []() consteval noexcept {
            if constexpr (Internal::isInvTrig<M1>())
                return Traits<M1>::Upper ? CblasUpper : CblasLower;
            else
                return Traits<M2>::Upper ? CblasUpper : CblasLower;
        }();

        constexpr auto Uplo = [](auto uploNoTrans) consteval noexcept {
            if constexpr (TransA == CblasTrans)
                return uploNoTrans == CblasUpper ? CblasLower : CblasUpper;
            return uploNoTrans;
        }(UploNoTrans);

        constexpr auto Diag = []() consteval noexcept {
            if constexpr (Internal::isInvTrig<M1>())
                return Traits<M1>::Unit ? CblasUnit : CblasNonUnit;
            else
                return Traits<M2>::Unit ? CblasUnit : CblasNonUnit;
        }();

        const size_t m = getRow();
        const size_t n = getCol();
        const auto* a = reinterpret_cast<const Tm*>([&]() noexcept {
            if constexpr (Internal::isInvTrig<M1>())
                return lhs.getExpr().getExpr().data_handle();
            else
                return rhs.getExpr().getExpr().data_handle();
        }());
        const size_t lda = [&]() noexcept {
            if constexpr (Internal::isInvTrig<M1>())
                return lhs.getExpr().getExpr().getMajorStride();
            else
                return rhs.getExpr().getExpr().getMajorStride();
        }();
        auto* const b = reinterpret_cast<Tm*>(target.data_handle());
        const size_t ldb = target.getMajorStride();
        if constexpr (isComplex()) {
            const Tc alpha = 1;
            if constexpr (T::Prec == Float32)
                cblas_ctrsm_64(Layout, Side, Uplo, TransA, Diag, m, n, (Tm*)&alpha, a, lda, b, ldb);
            else
                cblas_ztrsm_64(Layout, Side, Uplo, TransA, Diag, m, n, (Tm*)&alpha, a, lda, b, ldb);
        }
        else {
            const Tm alpha = 1;
            if constexpr (T::Prec == Float32)
                cblas_strsm_64(Layout, Side, Uplo, TransA, Diag, m, n, alpha, a, lda, b, ldb);
            else
                cblas_dtrsm_64(Layout, Side, Uplo, TransA, Diag, m, n, alpha, a, lda, b, ldb);
        }
    }
}
