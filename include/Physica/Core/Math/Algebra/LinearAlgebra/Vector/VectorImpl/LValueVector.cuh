/*
 * Copyright 2022-2023 WeiBo He.
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

#include "LValueVector.h"
#include "RValueVector.cuh"

namespace Physica::Core {
    template<class Derived>
    class device_obj<LValueVector<Derived>> : public device_obj<RValueVector<Derived>> {
    public:
        using Base = device_obj<RValueVector<Derived>>;
        using typename Base::ScalarType;
    public:
        ~device_obj() = default;
        /* Operators */
        __host__ __device__ inline device_obj& operator=(const device_obj& obj);
        __host__ __device__ inline device_obj& operator=(device_obj&& obj) noexcept;
        template<class OtherDerived>
        __host__ __device__ device_obj<Derived>& operator=(const device_obj<RValueVector<OtherDerived>>& v);
        template<class AnyScalar>
        inline device_obj<Derived>& operator=(const ScalarBase<AnyScalar>& s);
        [[nodiscard]] __device__ ScalarType& operator[](size_t index) { return *data_ptr(index); }
        [[nodiscard]] __device__ const ScalarType& operator[](size_t index) const { return *data_ptr(index); }
        /* Operations */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const { return *data_ptr(index); }
        /* Getters */
        [[nodiscard]] __host__ __device__ inline ScalarType* data_ptr(size_t index);
        [[nodiscard]] __host__ __device__ inline const ScalarType* data_ptr(size_t index) const;
    protected:
        device_obj() = default;
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
    };

    template<class Derived, class OtherDerived>
    __host__ __device__ inline void operator+=(device_obj<LValueVector<Derived>>& v1, const device_obj<RValueVector<OtherDerived>>& v2) {
        v1.getDerived() = v1.getDerived() + v2.getDerived();
    }

    template<class Derived, class OtherDerived>
    __host__ __device__ inline void operator-=(device_obj<LValueVector<Derived>>& v1, const device_obj<RValueVector<OtherDerived>>& v2) {
        v1.getDerived() = v1.getDerived() - v2.getDerived();
    }

    template<class VectorType, class ScalarType>
    __host__ __device__ inline void operator+=(device_obj<LValueVector<VectorType>>& v, const ScalarBase<ScalarType>& s) {
        v.getDerived() = v.getDerived() + s.getDerived();
    }

    template<class VectorType, class ScalarType>
    __host__ __device__ inline void operator*=(device_obj<LValueVector<VectorType>>& v, const ScalarBase<ScalarType>& s) {
        v.getDerived() = v.getDerived() * s.getDerived();
    }

    template<class VectorType, class ScalarType>
    __host__ __device__ inline void operator/=(device_obj<LValueVector<VectorType>>& v, const ScalarBase<ScalarType>& s) {
        v.getDerived() = v.getDerived() * reciprocal(s.getDerived());
    }
}

#include "LValueVectorImpl/LValueVectorImpl.cuh"
