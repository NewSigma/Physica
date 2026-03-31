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

#include "InvGEMM.h"
#include "Physica/Core/Exception/MKL/Lapack.h"

namespace Physica {
    template<Matrix M, bool Upper, bool Unit>
    void MatrixTrig<M, Upper, Unit>::assign_mkl(Matrix auto&& target) const {
        using Tm = decltype(std::declval<T>().toMKL());
        target.assert_assign_mkl(mat);
        target.zeros();

        constexpr auto Layout = MatrixMajor::isRowMatrix<M>() ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR;
        constexpr char Uplo = Upper ? 'U' : 'L';
        const size_t m = getRow();
        const size_t n = getCol();
        const auto* a = reinterpret_cast<const Tm*>(mat.data());
        const size_t lda = MatrixMajor::isRowMatrix<M>() ? n : m;
        auto* b = reinterpret_cast<Tm*>(target.data());
        const size_t ldb = lda;
        if constexpr (Base::isComplex()) {
            if constexpr (T::Prec == Float32)
                check_lapack(LAPACKE_clacpy_64(Layout, Uplo, m, n, a, lda, b, ldb));
            else
                check_lapack(LAPACKE_zlacpy_64(Layout, Uplo, m, n, a, lda, b, ldb));
        }
        else {
            if constexpr (T::Prec == Float32)
                check_lapack(LAPACKE_slacpy_64(Layout, Uplo, m, n, a, lda, b, ldb));
            else
                check_lapack(LAPACKE_dlacpy_64(Layout, Uplo, m, n, a, lda, b, ldb));
        }

        if constexpr (Unit)
            target.diag() = Trv(1);
    }
}
