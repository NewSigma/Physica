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

#include <mkl_cblas.h>
#include "../StridedVector.h"

namespace Physica {
    template<class Derived>
    void StridedVector<Derived>::assign_mkl(Vector auto& v) const noexcept {
        v.assert_assign_lapack(Base::getDerived());
        const size_t n = Base::getLength();
        const void* x = data_handle();
        const int strideX = static_cast<int>(getStride());
        const int strideY = static_cast<int>(v.getStride());
        void* y = v.data();
        if constexpr (T::isComplex()) {
            if constexpr (T::Prec == Float32)
                cblas_ccopy_64(n, x, strideX, y, strideY);
            else
                cblas_zcopy_64(n, x, strideX, y, strideY);
        }
        else {
            if constexpr (T::Prec == Float32)
                cblas_scopy_64(n, static_cast<const float*>(x), strideX, static_cast<float*>(y), strideY);
            else
                cblas_dcopy_64(n, static_cast<const double*>(x), strideX, static_cast<double*>(y), strideY);
        }
    }

    template<class Derived>
    auto StridedVector<Derived>::norm1_mkl() const noexcept -> Tr {
        const size_t n = Base::getLength();
        const void* x = data_handle();
        const int stride = static_cast<int>(getStride());
        if constexpr (T::isComplex()) {
            if constexpr (T::Prec == Float32)
                return cblas_scasum_64(n, x, stride);
            else {
                static_assert(T::Prec == Float64, "[Error]: Unexpected type");
                return cblas_dzasum_64(n, x, stride);
            }
        }
        else {
            if constexpr (T::Prec == Float32)
                return cblas_sasum_64(n, static_cast<const float*>(x), stride);
            else {
                static_assert(T::Prec == Float64, "[Error]: Unexpected type");
                return cblas_dasum_64(n, static_cast<const double*>(x), stride);
            }
        }
    }

    template<class Derived>
    auto StridedVector<Derived>::norm2_mkl() const noexcept -> Tr {
        const size_t n = Base::getLength();
        const void* x = data_handle();
        const int stride = static_cast<int>(getStride());
        if constexpr (T::isComplex()) {
            if constexpr (T::Prec == Float32)
                return cblas_scnrm2_64(n, x, stride);
            else {
                static_assert(T::Prec == Float64, "[Error]: Unexpected type");
                return cblas_dznrm2_64(n, x, stride);
            }
        }
        else {
            if constexpr (T::Prec == Float32)
                return cblas_snrm2_64(n, static_cast<const float*>(x), stride);
            else {
                static_assert(T::Prec == Float64, "[Error]: Unexpected type");
                return cblas_dnrm2_64(n, static_cast<const double*>(x), stride);
            }
        }
    }
}
