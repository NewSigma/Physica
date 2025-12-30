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

#include "Inverse.h"

namespace Physica {
    template<Scalar T, bool Pivot>
    void Inverse<DenseLU<T, Pivot>>::assign_mkl(Matrix auto& target) const {
        constexpr static int Layout = LAPACK_COL_MAJOR;
        target.assert_assign(*this);
        lu.getMatrixLU().assign(target);

        const size_t n = getRow();
        auto* a = reinterpret_cast<Tm*>(target.data());
        if constexpr (Pivot) {
            const MKL_INT64* ipiv = lu.getPerm().getIndices().data();
            if constexpr (isComplex) {
                if constexpr (T::Prec == Float32)
                    check_lapack(LAPACKE_cgetri_64(Layout, n, a, n, ipiv));
                else
                    check_lapack(LAPACKE_zgetri_64(Layout, n, a, n, ipiv));
            }
            else {
                if constexpr (T::Prec == Float32)
                    check_lapack(LAPACKE_sgetri_64(Layout, n, a, n, ipiv));
                else
                    check_lapack(LAPACKE_dgetri_64(Layout, n, a, n, ipiv));
            }
        }
        else {
            if constexpr (isComplex) {
                if constexpr (T::Prec == Float32)
                    check_lapack(LAPACKE_mkl_cgetrinp_64(Layout, n, a, n));
                else
                    check_lapack(LAPACKE_mkl_zgetrinp_64(Layout, n, a, n));
            }
            else {
                if constexpr (T::Prec == Float32)
                    check_lapack(LAPACKE_mkl_sgetrinp_64(Layout, n, a, n));
                else
                    check_lapack(LAPACKE_mkl_dgetrinp_64(Layout, n, a, n));
            }
        }
    }
}
