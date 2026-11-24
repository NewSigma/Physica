/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
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
        if constexpr (T::isComplex()) {
            if constexpr (T::Prec == Float32)
                cblas_cgemm_batch_64(CblasColMajor, &transA, &transB, &m, &n, &k, &alpha, (const void**)ptrs.data_ptr(0, 0), &lda, (const void**)ptrs.data_ptr(1, 0), &ldb, &beta, (void**)ptrs.data_ptr(2, 0), &ldc, groupCount, &groupSize);
            else {
                static_assert(T::Prec == Float64, "[Error]: Unknown ScalarType");
                cblas_zgemm_batch_64(CblasColMajor, &transA, &transB, &m, &n, &k, &alpha, (const void**)ptrs.data_ptr(0, 0), &lda, (const void**)ptrs.data_ptr(1, 0), &ldb, &beta, (void**)ptrs.data_ptr(2, 0), &ldc, groupCount, &groupSize);
            }
        }
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
