/*
 * Copyright 2023-2025 Weibo He.
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

#include "mkl_lapacke.h"
#include "../ContinuousMatrix.h"

namespace Physica {
    template<class Derived>
    auto ContinuousMatrix<Derived>::balance_mkl() -> VectorND<T> {
        Base::assert_balance();

        constexpr int Layout = MatrixMajor::getMajor<Derived>() == MatrixMajor::Row ? CblasRowMajor : CblasColMajor;
        constexpr char job = 'S';
        size_t n = Base::getRow();
        auto* a = reinterpret_cast<Tm*>(data());
        size_t lda = n;
        lapack_int ilo{}, ihi{};
        VectorND<T> result(n);
        auto* scale = reinterpret_cast<Tm*>(result.data());
        if constexpr (isComplex) {
            if constexpr (T::Prec == Float32)
                LAPACKE_cgebal(Layout, job, n, a, lda, &ilo, &ihi, scale);
            else
                LAPACKE_zgebal(Layout, job, n, a, lda, &ilo, &ihi, scale);
        }
        else {
            if constexpr (T::Prec == Float32)
                LAPACKE_sgebal(Layout, job, n, a, lda, &ilo, &ihi, scale);
            else
                LAPACKE_dgebal(Layout, job, n, a, lda, &ilo, &ihi, scale);
        }
        return reciprocal(result);
    }
}
