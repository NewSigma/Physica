/*
 * Copyright 2025-2026 Weibo He.
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
    template<Matrix M1, Matrix M2>
    void GEMM<M1, M2>::assign_mkl(Matrix auto&& target) const noexcept {
        using M = decltype(target);
        using Tm = decltype(std::declval<T>().toMKL());
        constexpr auto Layout = MatrixMajor::getMajor<M>() == MatrixMajor::Row ? CblasRowMajor : CblasColMajor;
        constexpr auto TransA = MatrixMajor::isSameMajor<M, M1>() ? CblasNoTrans : CblasTrans;
        constexpr auto TransB = MatrixMajor::isSameMajor<M, M2>() ? CblasNoTrans : CblasTrans;
        auto getData = [](const auto& mat) static {
            if constexpr (instanceof<Transpose, decltype(mat)>)
                return mat.transpose().data();
            else
                return mat.data();
        };

        const size_t m = getRow();
        const size_t n = getCol();
        const size_t k = mat1.getCol();
        const size_t lda = mat1.getMaxMinor();
        const size_t ldb = mat2.getMaxMinor();
        const size_t ldc = target.getMaxMinor();
        const auto* a = reinterpret_cast<const Tm*>(getData(mat1));
        const auto* b = reinterpret_cast<const Tm*>(getData(mat2));
        auto* c = reinterpret_cast<Tm*>(target.data());
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
