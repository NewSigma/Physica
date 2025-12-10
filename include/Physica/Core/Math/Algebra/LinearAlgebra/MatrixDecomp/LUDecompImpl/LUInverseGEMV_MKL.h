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

#include "LUInverseGEMV.h"

namespace Physica {
    template<Matrix M, Vector V> requires(instanceof_tx<Inverse, M> && requires { std::declval<M>().getLUDecomp(); })
    void GEMV<M, V>::assign_mkl(Vector auto& target) const {
        constexpr bool Pivot = Traits<M>::Pivot;
        constexpr int Layout = LAPACK_COL_MAJOR;
        constexpr char trans = 'N';
        size_t n = getLength();
        const auto* a = reinterpret_cast<const Tm*>(m.getLUDecomp().getMatrixLU().data());
        auto* b = reinterpret_cast<Tm*>(target.data());
        if constexpr (Pivot) {
            MKL_INT64* ipiv = m.getLUDecomp().getPerm().getIndices().data();
            if constexpr (isComplex) {
                if constexpr (T::Prec == Float32)
                    check_lapack(LAPACKE_cgetrs_work_64(Layout, trans, n, 1, a, n, ipiv, b, n));
                else
                    check_lapack(LAPACKE_zgetrs_work_64(Layout, trans, n, 1, a, n, ipiv, b, n));
            }
            else {
                if constexpr (T::Prec == Float32)
                    check_lapack(LAPACKE_sgetrs_work_64(Layout, trans, n, 1, a, n, ipiv, b, n));
                else
                    check_lapack(LAPACKE_dgetrs_work_64(Layout, trans, n, 1, a, n, ipiv, b, n));
            }
        }
        else {
            Array<MKL_INT64> perm(n);
            for (MKL_INT64 i = 0; i < n; ++i)
                perm[i] = i + 1;
            auto* ipiv = perm.data();
            if constexpr (isComplex) {
                if constexpr (T::Prec == Float32)
                    check_lapack(LAPACKE_cgetrs_work_64(Layout, trans, n, 1, a, n, ipiv, b, n));
                else
                    check_lapack(LAPACKE_zgetrs_work_64(Layout, trans, n, 1, a, n, ipiv, b, n));
            }
            else {
                if constexpr (T::Prec == Float32)
                    check_lapack(LAPACKE_sgetrs_work_64(Layout, trans, n, 1, a, n, ipiv, b, n));
                else
                    check_lapack(LAPACKE_dgetrs_work_64(Layout, trans, n, 1, a, n, ipiv, b, n));
            }
        }
    }
}
