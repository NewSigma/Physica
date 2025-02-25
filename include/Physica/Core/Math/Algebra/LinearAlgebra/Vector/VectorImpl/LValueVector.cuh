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

namespace Physica {
    template<class Derived>
    class device_obj<LValueVector<Derived>> : public device_obj<RValueVector<Derived>> {
        using This = device_obj<LValueVector<Derived>>;
        using Base = device_obj<RValueVector<Derived>>;
        template<size_t Length>
        using BlockType = device_obj<LVectorBlock<Derived, Length>>;
    public:
        using typename Base::ScalarType;
    protected:
        using PtrTy = ScalarType::PtrTy;
        using ConstPtrTy = ScalarType::ConstPtrTy;
        using RefTy = ScalarType::RefTy;
        using ConstRefTy = ScalarType::ConstRefTy;
    public:
        ~device_obj() = default;
        /* Operators */
        __host__ __device__ inline This& operator=(const This& obj);
        __host__ __device__ inline This& operator=(This&& obj);
        template<Vector V>
        __host__ __device__ device_obj<Derived>& operator=(const device_obj<V>& v);

        template<Scalar T> inline device_obj<Derived>& operator=(const T& x);
        template<Scalar T> __device__ void operator+=(const T& x);
        template<Scalar T> __device__ void operator-=(const T& x);
        template<Scalar T> __device__ void operator*=(const T& x);
        template<Scalar T> __device__ void operator/=(const T& x);

        template<Vector V> __device__ inline void operator+=(const device_obj<V>& v);
        template<Vector V> __device__ inline void operator-=(const device_obj<V>& v);
        [[nodiscard]] __device__ RefTy operator[](size_t index) { return *data_ptr(index); }
        [[nodiscard]] __device__ ConstRefTy operator[](size_t index) const { return *data_ptr(index); }
        /* Operations */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const { return *data_ptr(index); }

        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ inline auto head(size_t to) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ inline const auto head(size_t to) const noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ inline auto tail(size_t from) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ inline const auto tail(size_t from) const noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ inline auto segment(size_t from, size_t to) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ inline const auto segment(size_t from, size_t to) const noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ inline PtrTy data_ptr(size_t index);
        [[nodiscard]] __host__ __device__ inline ConstPtrTy data_ptr(size_t index) const;
    protected:
        device_obj() = default;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
    };
}

#include "LValueVectorImpl/LValueVectorImpl.cuh"
#include "LValueVectorImpl/LVectorBlock.h"
