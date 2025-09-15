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
    void VectorExpr<ExprType::Mul, V, U>::assign_mkl(Vector auto& v) const noexcept {
        using V1 = std::remove_cvref_t<decltype(v)>;
        assert(Base::getLength() == v.getLength());

        Base::getLHS().assign_check_mkl(v);
        if constexpr (std::same_as<std::remove_cvref_t<V>, V1>) {
            if (&Base::getLHS() != &v)
                getLHS().assign(v);
        }
        else
            getLHS().assign(v);

        const size_t n = Base::getLength();
        if constexpr (isComplex) {
            if constexpr (U::isComplex) {
                const auto a = T(Base::getRHS()).toMachine();
                if constexpr (T::Prec == Float32)
                    cblas_cscal_64(n, &a, v.data(), 1);
                else
                    cblas_zscal_64(n, &a, v.data(), 1);
            }
            else {
                const auto a = Base::getRHS().toMachine();
                if constexpr (T::Prec == Float32)
                    cblas_csscal_64(n, a, v.data(), 1);
                else
                    cblas_zdscal_64(n, a, v.data(), 1);
            }
        }
        else {
            const auto a = Base::getRHS().toMachine();
            if constexpr (T::Prec == Float32)
                cblas_sscal_64(n, a, reinterpret_cast<float*>(v.data()), 1);
            else
                cblas_dscal_64(n, a, reinterpret_cast<double*>(v.data()), 1);
        }
    }

    template<Vector V, Scalar U>
    void VectorExpr<ExprType::Mul, V, U>::assign_add_mkl(Vector auto& v) const noexcept {
        static_assert(T::Prec == Float32 || T::Prec == Float64);
        assert(Base::getLength() == v.getLength());
        const size_t n = Base::getLength();
        const auto alpha = Base::getRHS().toMachine();
        const void* x = Base::getLHS().data();
        if constexpr (isComplex) {
            if constexpr (T::Prec == Float32)
                cblas_caxpy_64(n, &alpha, x, 1, v.data(), 1);
            else
                cblas_zaxpy_64(n, &alpha, x, 1, v.data(), 1);
        }
        else {
            if constexpr (T::Prec == Float32)
                cblas_saxpy_64(n, alpha, static_cast<const float*>(x), 1, reinterpret_cast<float*>(v.data()), 1);
            else
                cblas_daxpy_64(n, alpha, static_cast<const double*>(x), 1, reinterpret_cast<double*>(v.data()), 1);
        }
    }
}
