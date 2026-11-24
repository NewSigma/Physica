/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "GEMV.h"

namespace Physica {
    template<Matrix M, Vector V>
    void GEMV<M, V>::assign_mkl(Vector auto&& target) const noexcept {
        using Tm = decltype(std::declval<T>().toMKL());
        static_assert(std::remove_cvref_t<decltype(mat)>::getMinorStrideAtCompile() == 1,
                      "[Error]: MKL expects a contiguous minor direction");
        constexpr auto Layout = mat.getMajor() == MatrixMajor::Row ? CblasRowMajor : CblasColMajor;
        constexpr auto Trans = CblasNoTrans;

        const size_t m = mat.getRow();
        const size_t n = mat.getCol();
        const size_t lda = [this]() {
            if constexpr (instanceof<M, Transpose>)
                return mat.transpose().getMajorStride();
            else
                return mat.getMajorStride();
        }();
        const auto* a = reinterpret_cast<const Tm*>([this]() {
            if constexpr (instanceof<M, Transpose>)
                return mat.transpose().data_handle();
            else
                return mat.data_handle();
        }());
        const auto* x = reinterpret_cast<const Tm*>(vec.data_handle());
        auto* y = reinterpret_cast<Tm*>(target.data_handle());
        const int incx = static_cast<int>(vec.getStride());
        const int incy = static_cast<int>(target.getStride());
        if constexpr (Base::isComplex()) {
            using Tc = T::ComplexType;
            const Tc alpha = 1;
            const Tc beta = 0;
            if constexpr (T::Prec == Float32)
                cblas_cgemv_64(Layout, Trans, m, n, (Tm*)&alpha, a, lda, x, incx, (Tm*)&beta, y, incy);
            else
                cblas_zgemv_64(Layout, Trans, m, n, (Tm*)&alpha, a, lda, x, incx, (Tm*)&beta, y, incy);
        }
        else {
            if constexpr (T::Prec == Float32)
                cblas_sgemv_64(Layout, Trans, m, n, 1, a, lda, x, incx, 0, y, incy);
            else
                cblas_dgemv_64(Layout, Trans, m, n, 1, a, lda, x, incx, 0, y, incy);
        }
    }
}
