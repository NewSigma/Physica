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

#include "GEMV.h"

namespace Physica {
    template<Matrix M, Vector V>
    void GEMV<M, V>::assign_mkl(Vector auto&& target) const noexcept {
        using Tm = decltype(std::declval<T>().toMKL());
        constexpr auto Layout = MatrixMajor::getMajor<M>() == MatrixMajor::Row ? CblasRowMajor : CblasColMajor;
        constexpr auto Trans = instanceof<Transpose, M> ? CblasTrans : CblasNoTrans;
        auto getData = [](const auto& mat) static {
            if constexpr (instanceof<Transpose, M>)
                return mat.transpose().data();
            else
                return mat.data();
        };

        const size_t m = mat.getRow();
        const size_t n = mat.getCol();
        const size_t lda = mat.getMaxMinor();
        const auto* a = reinterpret_cast<const Tm*>(getData(mat));
        const auto* x = reinterpret_cast<const Tm*>(vec.data());
        auto* y = reinterpret_cast<Tm*>(target.data());
        if constexpr (Base::isComplex()) {
            using Tc = T::ComplexType;
            const Tc alpha = 1;
            const Tc beta = 0;
            if constexpr (T::Prec == Float32)
                cblas_cgemv_64(Layout, Trans, m, n, (Tm*)&alpha, a, lda, x, 1, (Tm*)&beta, y, 1);
            else
                cblas_zgemv_64(Layout, Trans, m, n, (Tm*)&alpha, a, lda, x, 1, (Tm*)&beta, y, 1);
        }
        else {
            if constexpr (T::Prec == Float32)
                cblas_sgemv_64(Layout, Trans, m, n, 1, a, lda, x, 1, 0, y, 1);
            else
                cblas_dgemv_64(Layout, Trans, m, n, 1, a, lda, x, 1, 0, y, 1);
        }
    }
}
