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

#include "../Mul.h"

namespace Physica {
    template<Vector V, Scalar U>
    void VectorExpr<ExprType::Mul, V, U>::assign_add_mkl(void* y) const noexcept {
        const size_t n = Base::getLength();
        const auto alpha = Base::getRHS().toMachine();
        const void* x = Base::getLHS().data();
        if constexpr (isComplex) {
            if constexpr (T::Option == Float32)
                cblas_caxpy(n, &alpha, x, 1, y, 1);
            else {
                static_assert(T::Option == Float64, "[Error]: Unexpected type");
                cblas_zaxpy(n, &alpha, x, 1, y, 1);
            }
        }
        else {
            if constexpr (T::Option == Float32)
                cblas_saxpy(n, alpha, (float*)x, 1, (float*)y, 1);
            else {
                static_assert(T::Option == Float64, "[Error]: Unexpected type");
                cblas_daxpy(n, alpha, (double*)x, 1, (double*)y, 1);
            }
        }
    }
}
