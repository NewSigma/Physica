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

#include "TrigGEMM.h"

namespace Physica {
    template<Matrix M1, Matrix M2>
    void GEMM<TrigUpper<M1>, M2>::assign_mkl(Matrix auto& target) const noexcept {
        using M = std::remove_cvref_t<decltype(target)>;
        using Tc = T::ComplexType;
        using Tm = std::conditional<isComplex, typename Tc::MKL_Complex, typename T::MachineType>::type;
        constexpr auto Layout = MatrixOption::isRowMatrix<M>() ? CblasRowMajor : CblasColMajor;
        constexpr auto Side = CblasLeft;
        constexpr auto Uplo = CblasUpper;
        constexpr auto TransA = CblasNoTrans;
        constexpr auto Diag = CblasNonUnit;
        const M buffer = mat1;
        target = mat2;

        const size_t m = getRow();
        const size_t n = getCol();
        const Tm alpha = 1;
        const auto* a = reinterpret_cast<const Tm*>(buffer.data());
        const size_t lda = Side == CblasLeft ? m : n;
        auto* b = reinterpret_cast<Tm*>(target.data());
        const size_t ldb = Layout == CblasColMajor ? m : n;
        if constexpr (isComplex) {
            if constexpr (T::Prec == Float32)
                cblas_ctrmm_64(Layout, Side, Uplo, TransA, Diag, m, n, &alpha, a, lda, b, ldb);
            else
                cblas_ztrmm_64(Layout, Side, Uplo, TransA, Diag, m, n, &alpha, a, lda, b, ldb);
        }
        else {
            if constexpr (T::Prec == Float32)
                cblas_strmm_64(Layout, Side, Uplo, TransA, Diag, m, n, alpha, a, lda, b, ldb);
            else
                cblas_dtrmm_64(Layout, Side, Uplo, TransA, Diag, m, n, alpha, a, lda, b, ldb);
        }
    }
}
