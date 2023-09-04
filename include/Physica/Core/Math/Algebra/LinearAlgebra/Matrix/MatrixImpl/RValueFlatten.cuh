/*
 * Copyright 2023 WeiBo He.
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

#include "RValueFlatten.h"

namespace Physica::Core {
    namespace Internal {
        template<class T> class Traits;

        template<class MatrixType>
        class Traits<device_obj<RValueFlatten<MatrixType>>> : public Traits<RValueFlatten<MatrixType>> {};
    }

    template<class MatrixType>
    class device_obj<RValueFlatten<MatrixType>> : public device_obj<RValueVector<RValueFlatten<MatrixType>>> {
        const device_obj<MatrixType>& mat;
    public:
        using host_obj = RValueFlatten<MatrixType>;
        using Base = device_obj<RValueVector<host_obj>>;
        using typename Base::ScalarType;
    public:
        __host__ __device__ device_obj(const device_obj<RValueMatrix<MatrixType>>& mat_) : mat(mat_.getDerived()) {}
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const;
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return mat.getRow() * mat.getColumn(); }
    };

    template<class MatrixType>
    __device__ typename device_obj<RValueFlatten<MatrixType>>::ScalarType device_obj<RValueFlatten<MatrixType>>::calc(size_t index) const {
        const size_t major = index / mat.getMaxMinor();
        const size_t minor = index % mat.getMaxMinor();
        return mat.calcFromMajorMinor(major, minor);
    }
}
