/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <mkl_trans.h>
#include "Transpose.h"

namespace Physica {
    template<Matrix M> requires(std::remove_cvref_t<M>::isCompact())
    void Transpose<M>::assign_mkl(Matrix auto&& target) const {
        static_assert(MatrixMajor::isSameMajor<M, decltype(target)>(), "[Error]: Cannot apply MKL to this expr");
        static_assert(target.getMinorStrideAtCompile() == 1, "[Error]: MKL expects a contiguous minor direction");
        using Tm = decltype(std::declval<T>().toMKL());
        target.assert_assign_lapack(mat);

        constexpr char ordering = MatrixMajor::isRowMatrix<M>() ? 'R' : 'C';
        constexpr char trans = 'T';
        const size_t rows = getRow();
        const size_t cols = getCol();
        const auto* A = reinterpret_cast<const Tm*>(mat.data_handle());
        auto* B = reinterpret_cast<Tm*>(target.data_handle());
        const size_t lda = mat.getMajorStride();
        const size_t ldb = target.getMajorStride();
        if constexpr (T::isComplex()) {
            const auto alpha = T(1).toMKL();
            if constexpr (T::Prec == Float32)
                mkl_comatcopy(ordering, trans, rows, cols, alpha, A, lda, B, ldb);
            else
                mkl_zomatcopy(ordering, trans, rows, cols, alpha, A, lda, B, ldb);
        }
        else {
            constexpr double alpha = 1;
            if constexpr (T::Prec == Float32)
                mkl_somatcopy(ordering, trans, rows, cols, alpha, A, lda, B, ldb);
            else
                mkl_domatcopy(ordering, trans, rows, cols, alpha, A, lda, B, ldb);
        }
    }
}
