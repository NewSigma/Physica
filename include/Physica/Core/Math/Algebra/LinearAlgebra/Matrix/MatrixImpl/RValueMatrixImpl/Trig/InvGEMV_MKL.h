/*
 * Copyright 2026 Weibo He.
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

#include "InvGEMV.h"
#include "Physica/Core/Exception/MKL/Lapack.h"

namespace Physica {
    template<Matrix M, Vector V> requires(instanceof<M, Inverse> && instanceof_tx<typename Traits<M>::ExprType, MatrixTrig>)
    void GEMV<M, V>::assign_mkl(Vector auto& target) const {
        using Tm = decltype(std::declval<T>().toMKL());
        rhs.assign(target);

        const auto& mat = inv.getExpr().getExpr();
        constexpr auto Layout = MatrixMajor::isRowMatrix<decltype(mat)>() ? CblasRowMajor : CblasColMajor;
        constexpr auto Uplo = Traits<M>::Upper ? CblasUpper : CblasLower;
        constexpr auto TransA = CblasNoTrans;
        constexpr auto Diag = Traits<M>::Unit ? CblasUnit : CblasNonUnit;

        const size_t n = getLength();
        const auto* a = reinterpret_cast<const Tm*>(mat.data_handle());
        const size_t lda = mat.getMajorStride();
        auto* x = reinterpret_cast<Tm*>(target.data_handle());
        if constexpr (isComplex()) {
            if constexpr (T::Prec == Float32)
                cblas_ctrsv_64(Layout, Uplo, TransA, Diag, n, reinterpret_cast<const void*>(a), lda, reinterpret_cast<void*>(x), 1);
            else
                cblas_ztrsv_64(Layout, Uplo, TransA, Diag, n, reinterpret_cast<const void*>(a), lda, reinterpret_cast<void*>(x), 1);
        }
        else {
            if constexpr (T::Prec == Float32)
                cblas_strsv_64(Layout, Uplo, TransA, Diag, n, a, lda, x, 1);
            else
                cblas_dtrsv_64(Layout, Uplo, TransA, Diag, n, a, lda, x, 1);
        }
    }
}
