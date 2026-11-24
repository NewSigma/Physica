/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "GEMM.h"

namespace Physica {
    template<Matrix M1, Matrix M2>
    void GEMM<M1, M2>::assign_mkl(Matrix auto&& target) const noexcept {
        using M = decltype(target);
        using Tm = decltype(std::declval<T>().toMKL());
        static_assert(isMinorContiguous(target), "[Error]: MKL expects a contiguous minor direction");
        constexpr auto Layout = target.getMajor() == MatrixMajor::Row ? CblasRowMajor : CblasColMajor;
        constexpr auto TransA = MatrixMajor::isSameMajor<M, M1>() ? CblasNoTrans : CblasTrans;
        constexpr auto TransB = MatrixMajor::isSameMajor<M, M2>() ? CblasNoTrans : CblasTrans;
        auto getData = [](const auto& mat) static {
            if constexpr (instanceof<decltype(mat), Transpose>)
                return mat.transpose().data_handle();
            else
                return mat.data_handle();
        };
        auto getMajorStride = [](const auto& mat) static {
            if constexpr (instanceof<decltype(mat), Transpose>)
                return mat.transpose().getMajorStride();
            else
                return mat.getMajorStride();
        };

        const size_t m = getRow();
        const size_t n = getCol();
        const size_t k = mat1.getCol();
        const size_t lda = getMajorStride(mat1);
        const size_t ldb = getMajorStride(mat2);
        const size_t ldc = target.getMajorStride();
        const auto* a = reinterpret_cast<const Tm*>(getData(mat1));
        const auto* b = reinterpret_cast<const Tm*>(getData(mat2));
        auto* c = reinterpret_cast<Tm*>(target.data_handle());
        if constexpr (Base::isComplex()) {
            using Tc = T::ComplexType;
            const Tc alpha = 1;
            const Tc beta = 0;
            if constexpr (T::Prec == Float32)
                cblas_cgemm_64(Layout, TransA, TransB, m, n, k, (Tm*)&alpha, a, lda, b, ldb, (Tm*)&beta, c, ldc);
            else
                cblas_zgemm_64(Layout, TransA, TransB, m, n, k, (Tm*)&alpha, a, lda, b, ldb, (Tm*)&beta, c, ldc);
        }
        else {
            if constexpr (T::Prec == Float32)
                cblas_sgemm_64(Layout, TransA, TransB, m, n, k, 1, a, lda, b, ldb, 0, c, ldc);
            else
                cblas_dgemm_64(Layout, TransA, TransB, m, n, k, 1, a, lda, b, ldb, 0, c, ldc);
        }
    }
}
