/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "InvGEMV.h"
#include "Physica/Core/Exception/MKL/Lapack.h"

namespace Physica {
    template<Matrix M, Vector V> requires(instanceof<M, Inverse> && instanceof_tx<typename Traits<M>::ExprType, MatrixTrig>)
    void GEMV<M, V>::assign_mkl(Vector auto& target) const {
        using Tm = decltype(std::declval<T>().toMKL());
        using Data = std::remove_cvref_t<decltype(std::declval<const std::remove_cvref_t<M>&>().getExpr().getExpr())>;
        static_assert(Data::getMinorStrideAtCompile() == 1, "[Error]: MKL expects a contiguous minor direction");
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
        const int incx = static_cast<int>(target.getStride());
        if constexpr (isComplex()) {
            if constexpr (T::Prec == Float32)
                cblas_ctrsv_64(Layout, Uplo, TransA, Diag, n, reinterpret_cast<const void*>(a), lda, reinterpret_cast<void*>(x), incx);
            else
                cblas_ztrsv_64(Layout, Uplo, TransA, Diag, n, reinterpret_cast<const void*>(a), lda, reinterpret_cast<void*>(x), incx);
        }
        else {
            if constexpr (T::Prec == Float32)
                cblas_strsv_64(Layout, Uplo, TransA, Diag, n, a, lda, x, incx);
            else
                cblas_dtrsv_64(Layout, Uplo, TransA, Diag, n, a, lda, x, incx);
        }
    }
}
