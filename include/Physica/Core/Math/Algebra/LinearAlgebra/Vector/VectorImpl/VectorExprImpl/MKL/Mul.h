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
        Base::getLHS().assert_assign_mkl(v);

        using V1 = std::remove_cvref_t<decltype(v)>;
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
    void VectorExpr<ExprType::Mul, V, U>::assign_add_mkl(Vector auto&& v) const noexcept {
        Base::getLHS().assert_assign_mkl(v);

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

    template<Vector V1, Vector V2>
    void VectorExpr<ExprType::Mul, V1, V2>::assign_mkl(Vector auto& v) const noexcept {
        Base::getLHS().assert_assign_mkl(v);
        Base::getRHS().assert_assign_mkl(v);

        size_t n = Base::getLength();
        const auto* a = reinterpret_cast<const Tm*>(Base::getLHS().data());
        const auto* b = reinterpret_cast<const Tm*>(Base::getRHS().data());
        auto* r = reinterpret_cast<Tm*>(v.data());
        if constexpr (isComplex) {
            if constexpr (T::Prec == Float32)
                vcMul_64(n, a, b, r);
            else
                vzMul_64(n, a, b, r);
        }
        else {
            if constexpr (T::Prec == Float16)
                vhMul_64(n, a, b, r);
            else if constexpr (T::Prec == Float32)
                vsMul_64(n, a, b, r);
            else
                vdMul_64(n, a, b, r);
        }
    }
}
