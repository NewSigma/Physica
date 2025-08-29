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

#include "Physica/Core/Exception/MKL/Lapack.h"
#include "LUDecomp.h"

namespace Physica {
    template<Scalar T, bool Pivot>
    void LUDecomp<T, Pivot>::compute_mkl(const Matrix auto& source) {
        pre_compute(source);
        working = source;

        constexpr int Layout = LAPACK_COL_MAJOR;
        const size_t m = source.getRow();
        const size_t n = source.getCol();
        const size_t lda = m;
        auto* a = reinterpret_cast<Tm*>(working.data());
        if constexpr (Pivot) {
            auto* ipiv = reinterpret_cast<MKL_INT64*>(biasOrder.data());
            if constexpr (isComplex) {
                if constexpr (T::Prec == Float32)
                    check_lapack(LAPACKE_cgetrf_64(Layout, m, n, a, lda, ipiv));
                else
                    check_lapack(LAPACKE_zgetrf_64(Layout, m, n, a, lda, ipiv));
            }
            else {
                if constexpr (T::Prec == Float32)
                    check_lapack(LAPACKE_sgetrf_64(Layout, m, n, a, lda, ipiv));
                else
                    check_lapack(LAPACKE_dgetrf_64(Layout, m, n, a, lda, ipiv));
            }
        }
        else {
            if constexpr (isComplex) {
                if constexpr (T::Prec == Float32)
                    check_lapack(LAPACKE_mkl_cgetrfnp_64(Layout, m, n, a, lda));
                else
                    check_lapack(LAPACKE_mkl_zgetrfnp_64(Layout, m, n, a, lda));
            }
            else {
                if constexpr (T::Prec == Float32)
                    check_lapack(LAPACKE_mkl_sgetrfnp_64(Layout, m, n, a, lda));
                else
                    check_lapack(LAPACKE_mkl_dgetrfnp_64(Layout, m, n, a, lda));
            }
        }
    }
}
