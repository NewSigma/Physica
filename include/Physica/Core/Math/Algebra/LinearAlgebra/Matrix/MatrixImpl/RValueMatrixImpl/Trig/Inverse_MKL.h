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

#include "Inverse.h"
#include "Physica/Core/Exception/MKL/Lapack.h"

namespace Physica {
    template<Matrix M> requires(instanceof_tx<M, MatrixTrig>)
    void Inverse<M>::assign_mkl(Matrix auto& target) const {
        using Tm = decltype(std::declval<T>().toMKL());
        trig.assign(target);

        constexpr auto Layout = MatrixMajor::isRowMatrix<decltype(target)>() ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR;
        constexpr char Uplo = Traits<M>::Upper ? 'U' : 'L';
        constexpr char Diag = Traits<M>::Unit ? 'U' : 'N';
        size_t n = getRow();
        auto* a = reinterpret_cast<Tm*>(target.data());
        if constexpr (Base::isComplex()) {
            if constexpr (T::Prec == Float32)
                check_lapack(LAPACKE_ctrtri_64(Layout, Uplo, Diag, n, a, n));
            else
                check_lapack(LAPACKE_ztrtri_64(Layout, Uplo, Diag, n, a, n));
        }
        else {
            if constexpr (T::Prec == Float32)
                check_lapack(LAPACKE_strtri_64(Layout, Uplo, Diag, n, a, n));
            else
                check_lapack(LAPACKE_dtrtri_64(Layout, Uplo, Diag, n, a, n));
        }
    }
}
