/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "SymmEigenSolver.h"
#include "Physica/Core/Exception/MKL/Lapack.h"

namespace Physica {
    template<Scalar T, size_t Order>
    void SymmEigenSolver<T, Order>::compute_mkl(const Matrix auto& source) {
        static_assert(T::Prec == Float32 || T::Prec == Float64);
        using Tm = decltype(std::declval<T>().toMKL());
        constexpr int Major = MatrixMajor::isColMatrix<decltype(source)>() ? MatrixMajor::Col : MatrixMajor::Row;
        constexpr int Layout = Major == MatrixMajor::Col ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;

        eigenvectors = source;

        const int jobz = getNeedEigenvectors() ? 'V' : 'N';
        const size_t n = source.getRow();
        const size_t lda = n;
        auto* a = reinterpret_cast<Tm*>(eigenvectors.data());
        auto* w = reinterpret_cast<Tm*>(eigenvalues.data());
        if constexpr (T::Prec == Float32)
            check_lapack(LAPACKE_ssyevd_64(Layout, jobz, 'U', n, a, lda, w));
        else
            check_lapack(LAPACKE_dsyevd_64(Layout, jobz, 'U', n, a, lda, w));
    }
}
