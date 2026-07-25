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
    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig> || instanceof_tx<M2, MatrixTrig>)
    void GEMM<M1, M2>::assign_mkl(Matrix auto& target) const noexcept {
        using M = std::remove_cvref_t<decltype(target)>;
        using Tm = decltype(std::declval<T>().toMKL());
        constexpr bool TrigLHS = instanceof_tx<M1, MatrixTrig>;
        using Trig = std::conditional_t<TrigLHS, M1, M2>;
        constexpr auto Layout = MatrixMajor::isRowMatrix<M>() ? CblasRowMajor : CblasColMajor;
        constexpr auto Side = TrigLHS ? CblasLeft : CblasRight;
        constexpr auto TransTrig = MatrixMajor::isSameMajor<Trig, M>() ? CblasNoTrans : CblasTrans;
        constexpr auto UploNoTrans = Traits<Trig>::Upper ? CblasUpper : CblasLower;
        constexpr auto Uplo = []() consteval noexcept {
            if constexpr (TransTrig == CblasTrans)
                return UploNoTrans == CblasUpper ? CblasLower : CblasUpper;
            return UploNoTrans;
        }();

        constexpr auto Diag = Traits<Trig>::Unit ? CblasUnit : CblasNonUnit;

        if constexpr (TrigLHS)
            rhs.assign(target);
        else
            lhs.assign(target);

        const size_t m = getRow();
        const size_t n = getCol();
        const auto* a = reinterpret_cast<const Tm*>([this]() noexcept {
            if constexpr (TrigLHS)
                return lhs.getExpr().data_handle();
            else
                return rhs.getExpr().data_handle();
        }());
        const size_t lda = [this]() noexcept {
            if constexpr (TrigLHS)
                return lhs.getExpr().getMajorStride();
            else
                return rhs.getExpr().getMajorStride();
        }();
        auto* b = reinterpret_cast<Tm*>(target.data_handle());
        const size_t ldb = target.getMajorStride();
        if constexpr (Base::isComplex()) {
            const Tc alpha = 1;
            if constexpr (T::Prec == Float32)
                cblas_ctrmm_64(Layout, Side, Uplo, TransTrig, Diag, m, n, (Tm*)&alpha, a, lda, b, ldb);
            else
                cblas_ztrmm_64(Layout, Side, Uplo, TransTrig, Diag, m, n, (Tm*)&alpha, a, lda, b, ldb);
        }
        else {
            const Tm alpha = 1;
            if constexpr (T::Prec == Float32)
                cblas_strmm_64(Layout, Side, Uplo, TransTrig, Diag, m, n, alpha, a, lda, b, ldb);
            else
                cblas_dtrmm_64(Layout, Side, Uplo, TransTrig, Diag, m, n, alpha, a, lda, b, ldb);
        }
    }
}
