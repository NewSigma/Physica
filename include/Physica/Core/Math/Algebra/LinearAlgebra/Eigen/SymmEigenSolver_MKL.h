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

#include "SymmEigenSolver.h"
#include "Physica/Core/Exception/MKL/Lapack.h"

namespace Physica {
    template<Scalar T, size_t Order>
    void SymmEigenSolver<T, Order>::compute_mkl(const Matrix auto& source) {
        static_assert(T::Prec == Float32 || T::Prec == Float64);
        constexpr int Major = MatrixOption::isColMatrix<decltype(source)>() ? MatrixOption::Col : MatrixOption::Row;
        constexpr int Layout = Major == MatrixOption::Col ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;

        eigenvectors = source;

        const int jobz = getNeedEigenvectors() ? 'V' : 'N';
        const size_t n = source.getRow();
        const size_t lda = n;
        auto* a = reinterpret_cast<Tm*>(eigenvectors.data());
        auto* w = reinterpret_cast<Tm*>(eigenvalues.data());
        if constexpr (T::Prec == Float32)
            check_lapack(LAPACKE_ssyev_64(Layout, jobz, 'U', n, a, lda, w));
        else
            check_lapack(LAPACKE_dsyev_64(Layout, jobz, 'U', n, a, lda, w));
    }
}
