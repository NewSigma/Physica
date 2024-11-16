/*
 * Copyright 2022-2024 Weibo He.
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
        using Base = device_obj<RValueVector<Derived>>;
        template<size_t Length>
        using BlockType = device_obj<LVectorBlock<Derived, Length>>;
    public:
        using typename Base::ScalarType;
    protected:
        using PtrTy = typename ScalarType::PtrTy;
        using ConstPtrTy = typename ScalarType::ConstPtrTy;
        using RefTy = typename ScalarType::RefTy;
        using ConstRefTy = typename ScalarType::ConstRefTy;
    public:
        ~device_obj() = default;
        /* Operators */
        __host__ __device__ inline device_obj& operator=(const device_obj& obj);
        __host__ __device__ inline device_obj& operator=(device_obj&& obj);
        template<Vector V>
        __host__ __device__ device_obj<Derived>& operator=(const device_obj<V>& v);

        template<Scalar T>
        inline device_obj<Derived>& operator=(const T& s);
        void operator+=(const ScalarType& s) { (*this) = (*this) + s; }
        void operator-=(const ScalarType& s) { (*this) = (*this) - s; }
        void operator*=(const ScalarType& s) { (*this) = (*this) * s; }
        void operator/=(const ScalarType& s) { (*this) = (*this) / s; }

        [[nodiscard]] __device__ RefTy operator[](size_t index) { return *data_ptr(index); }
        [[nodiscard]] __device__ ConstRefTy operator[](size_t index) const { return *data_ptr(index); }
        /* Operations */
        template<Side Owner = GetSide()>
        [[nodiscard]] __device__ ScalarType calc(size_t index) const { return *data_ptr(index); }

        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ inline BlockType<Length> head(size_t to);
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ inline const BlockType<Length> head(size_t to) const;
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ inline BlockType<Length> tail(size_t from);
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ inline const BlockType<Length> tail(size_t from) const;
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ inline BlockType<Length> segment(size_t from, size_t to);
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ inline const BlockType<Length> segment(size_t from, size_t to) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ inline PtrTy data_ptr(size_t index);
        [[nodiscard]] __host__ __device__ inline ConstPtrTy data_ptr(size_t index) const;
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

    template<class VectorType, Scalar T>
    __host__ __device__ inline void operator+=(device_obj<LValueVector<VectorType>>& v, const T& x) {
        v.getDerived() = v.getDerived() + x;
    }

    template<class VectorType, Scalar T>
    __host__ __device__ inline void operator*=(device_obj<LValueVector<VectorType>>& v, const T& x) {
        v.getDerived() = v.getDerived() * x;
    }

    template<class VectorType, Scalar T>
    __host__ __device__ inline void operator/=(device_obj<LValueVector<VectorType>>& v, const T& x) {
        v.getDerived() = v.getDerived() * reciprocal(x);
    }
}

#include "LValueVectorImpl/LValueVectorImpl.cuh"
#include "LValueVectorImpl/LVectorBlock.h"
