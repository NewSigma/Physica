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

#include <mkl_trans.h>
#include "Transpose.h"

namespace Physica {
    template<Matrix M> requires(std::remove_cvref_t<M>::isCompact())
    void Transpose<M>::assign_mkl(Matrix auto&& target) const {
        static_assert(MatrixMajor::isSameMajor<M, decltype(target)>(), "[Error]: Cannot apply MKL to this expr");
        using Tm = decltype(std::declval<T>().toMKL());
        target.assert_assign_lapack(mat);
        
        constexpr char ordering = MatrixMajor::isRowMatrix<M>() ? 'R' : 'C';
        constexpr char trans = 'T';
        const size_t rows = getRow();
        const size_t cols = getCol();
        const auto* A = reinterpret_cast<const Tm*>(mat.data());
        auto* B = reinterpret_cast<Tm*>(target.data());
        const size_t lda = Base::getMaxMinor();
        const size_t ldb = target.getMaxMinor();
        if constexpr (T::isComplex()) {
            const auto alpha = T(1).toMKL();
            if constexpr (T::Prec == Float32)
                mkl_comatcopy(ordering, trans, rows, cols, alpha, A, lda, B, ldb);
            else
                mkl_zomatcopy(ordering, trans, rows, cols, alpha, A, lda, B, ldb);
        }
        else {
            constexpr double alpha = 1;
            if constexpr (T::Prec == Float32)
                mkl_somatcopy(ordering, trans, rows, cols, alpha, A, lda, B, ldb);
            else
                mkl_domatcopy(ordering, trans, rows, cols, alpha, A, lda, B, ldb);
        }
    }
}
