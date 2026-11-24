/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "InvGEMM.h"
#include "Physica/Core/Exception/MKL/Lapack.h"

namespace Physica {
    template<Matrix M, bool Upper, bool Unit>
    void MatrixTrig<M, Upper, Unit>::assign_mkl(Matrix auto&& target) const {
        static_assert(mat.getMinorStrideAtCompile() == 1 && target.getMinorStrideAtCompile() == 1, "[Error]: LAPACK expects a contiguous minor direction");
        using Tm = decltype(std::declval<T>().toMKL());
        target.assert_assign_lapack(mat);
        target.zeros();

        constexpr auto Layout = MatrixMajor::isRowMatrix<M>() ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR;
        constexpr char Uplo = Upper ? 'U' : 'L';
        const size_t m = getRow();
        const size_t n = getCol();
        const auto* a = reinterpret_cast<const Tm*>(mat.data_handle());
        const size_t lda = mat.getMajorStride();
        auto* b = reinterpret_cast<Tm*>(target.data_handle());
        const size_t ldb = target.getMajorStride();
        if constexpr (Base::isComplex()) {
            if constexpr (T::Prec == Float32)
                check_lapack(LAPACKE_clacpy_64(Layout, Uplo, m, n, a, lda, b, ldb));
            else
                check_lapack(LAPACKE_zlacpy_64(Layout, Uplo, m, n, a, lda, b, ldb));
        }
        else {
            if constexpr (T::Prec == Float32)
                check_lapack(LAPACKE_slacpy_64(Layout, Uplo, m, n, a, lda, b, ldb));
            else
                check_lapack(LAPACKE_dlacpy_64(Layout, Uplo, m, n, a, lda, b, ldb));
        }

        if constexpr (Unit)
            target.diag() = Trv(1);
    }
}
