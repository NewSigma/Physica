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

#include <mkl_cblas.h>
#include "BatchedGEMM.h"

namespace Physica {
    template<Scalar T>
    void BatchedGEMM<T>::compute_mkl(T alpha_, T beta_) noexcept {
        using Tm = decltype(std::declval<T>().toMKL());
        const auto transA = shape.transA ? CblasTrans : CblasNoTrans;
        const auto transB = shape.transB ? CblasTrans : CblasNoTrans;
        const auto m = static_cast<MKL_INT64>(shape.row);
        const auto n = static_cast<MKL_INT64>(shape.col);
        const auto k = static_cast<MKL_INT64>(shape.inner);
        const auto lda = static_cast<MKL_INT64>(shape.majorStrideA);
        const auto ldb = static_cast<MKL_INT64>(shape.majorStrideB);
        const auto ldc = static_cast<MKL_INT64>(shape.majorStrideC);
        const Tm alpha = alpha_.toMKL();
        const Tm beta = beta_.toMKL();
        const auto groupCount = static_cast<MKL_INT64>(1);
        const auto groupSize = static_cast<MKL_INT64>(getBatchSize());
        if constexpr (T::isComplex())
            cblas_cgemm_batch_64(CblasColMajor, &transA, &transB, &m, &n, &k, &alpha, (const void**)ptrs.data_ptr(0, 0), &lda, (const void**)ptrs.data_ptr(1, 0), &ldb, &beta, (void**)ptrs.data_ptr(2, 0), &ldc, groupCount, &groupSize);
        else {
            if constexpr (T::Prec == Float32)
                cblas_sgemm_batch_64(CblasColMajor, &transA, &transB, &m, &n, &k, &alpha, (const Tm**)ptrs.data_ptr(0, 0), &lda, (const Tm**)ptrs.data_ptr(1, 0), &ldb, &beta, (Tm**)ptrs.data_ptr(2, 0), &ldc, groupCount, &groupSize);
            else {
                static_assert(T::Prec == Float64, "[Error]: Unknown ScalarType");
                cblas_dgemm_batch_64(CblasColMajor, &transA, &transB, &m, &n, &k, &alpha, (const Tm**)ptrs.data_ptr(0, 0), &lda, (const Tm**)ptrs.data_ptr(1, 0), &ldb, &beta, (Tm**)ptrs.data_ptr(2, 0), &ldc, groupCount, &groupSize);
            }
        }
    }
}
