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

#include "GEMM.h"

namespace Physica {
    template<Matrix T1, Matrix T2>
    template<Matrix M>
    void MatrixProduct<T1, T2>::assign_mkl(LValueMatrix<M>& target) const {
        constexpr int Major = MatrixOption::isRowMatrix<T1>() ? MatrixOption::Row : MatrixOption::Col;
        constexpr auto Layout = Major == MatrixOption::Row ? CblasRowMajor : CblasColMajor;
        using ComplexType = ScalarType::ComplexType;
        using Tm = std::conditional<isComplex, typename ComplexType::MKL_Complex, typename ScalarType::MachineType>::type;
        const size_t m = getRow();
        const size_t n = getCol();
        const size_t k = mat1.getCol();
        const size_t lda = mat1.getMaxMinor();
        const size_t ldb = mat2.getMaxMinor();
        const size_t ldc = target.getMaxMinor();
        auto* a = reinterpret_cast<const Tm*>(mat1.data());
        auto* b = reinterpret_cast<const Tm*>(mat2.data());
        auto* c = reinterpret_cast<Tm*>(target.getDerived().data());

        if constexpr (isComplex) {
            const Tm alpha = 1;
            const Tm beta = 0;
            if constexpr (ScalarType::Option == Float32)
                cblas_cgemm_64(Layout, CblasNoTrans, CblasNoTrans, m, n, k, &alpha, a, lda, b, ldb, &beta, c, ldc);
            else
                cblas_zgemm_64(Layout, CblasNoTrans, CblasNoTrans, m, n, k, &alpha, a, lda, b, ldb, &beta, c, ldc);
        }
        else {
            if constexpr (ScalarType::Option == Float32)
                cblas_sgemm_64(Layout, CblasNoTrans, CblasNoTrans, m, n, k, 1, a, lda, b, ldb, 0, c, ldc);
            else
                cblas_dgemm_64(Layout, CblasNoTrans, CblasNoTrans, m, n, k, 1, a, lda, b, ldb, 0, c, ldc);
        }
    }
}
