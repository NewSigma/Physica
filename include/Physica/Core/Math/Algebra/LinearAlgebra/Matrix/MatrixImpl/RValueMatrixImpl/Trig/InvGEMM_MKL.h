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
    template<Matrix M1, Matrix M2> requires(instanceof<Inverse, M1> && instanceof_tx<MatrixTrig, typename Traits<M1>::ExprType>)
    void GEMM<M1, M2>::assign_mkl(Matrix auto& target) const noexcept {
        using M = std::remove_cvref_t<decltype(target)>;
        constexpr auto Layout = MatrixOption::isRowMatrix<M>() ? CblasRowMajor : CblasColMajor;
        constexpr auto Side = CblasLeft;
        constexpr auto Uplo = Traits<M1>::Upper ? CblasUpper : CblasLower;
        constexpr auto TransA = CblasNoTrans;
        constexpr auto Diag = Traits<M1>::Unit ? CblasUnit : CblasNonUnit;
        const M buffer = inv.getExpr();
        rhs.assign(target);

        const size_t m = getRow();
        const size_t n = getCol();
        const auto* a = reinterpret_cast<const Tm*>(buffer.data());
        const size_t lda = Side == CblasLeft ? m : n;
        auto* b = reinterpret_cast<Tm*>(target.data());
        const size_t ldb = Layout == CblasColMajor ? m : n;
        if constexpr (isComplex) {
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
