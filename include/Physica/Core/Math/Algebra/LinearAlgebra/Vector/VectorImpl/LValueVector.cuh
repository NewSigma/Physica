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
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
        using PtrTy = T::PtrTy;
        using ConstPtrTy = T::ConstPtrTy;
        using RefTy = T::RefTy;
        using ConstRefTy = T::ConstRefTy;
    public:
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This& obj) = delete;
        This& operator=(This&& obj) = delete;

        template<Scalar U>
        device_obj<Derived>& operator=(const U& x);
        __host__ __device__ void operator+=(const Scalar auto& x);
        __host__ __device__ void operator-=(const Scalar auto& x);
        __host__ __device__ void operator*=(const Scalar auto& x);
        __host__ __device__ void operator/=(const Scalar auto& x);

        __host__ __device__ device_obj<Derived>& operator=(const Vector auto& v);
        __host__ __device__ void operator+=(const Vector auto& v);
        __host__ __device__ void operator-=(const Vector auto& v);

        [[nodiscard]] __device__ RefTy operator[](size_t index);
        [[nodiscard]] __device__ ConstRefTy operator[](size_t index) const;
        /* Operations */
        [[nodiscard]] __device__ ConstRefTy calc(size_t index) const { return operator[](index); }
        [[nodiscard]] __device__ Tv calc_value(size_t index) const { return calc(index).value(); }

        __host__ __device__ void reverse(const auto& grad) const noexcept;

        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ auto head(size_t to) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ const auto head(size_t to) const noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ auto tail(size_t from) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ const auto tail(size_t from) const noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ auto segment(size_t from, size_t to) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ const auto segment(size_t from, size_t to) const noexcept;

        template<Matrix M>
        [[nodiscard]] __host__ __device__ auto reshape(const M& mat) noexcept;
        template<Matrix M>
        [[nodiscard]] __host__ __device__ const auto reshape(const M& mat) const noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ auto reshape_col(size_t row, size_t col) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ const auto reshape_col(size_t row, size_t col) const noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ auto reshape_row(size_t row, size_t col) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ const auto reshape_row(size_t row, size_t col) const noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ PtrTy data_ptr(size_t index) noexcept;
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr(size_t index) const noexcept;
    protected:
        device_obj() = default;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
    };
}

#include "LValueVectorImpl/LValueVectorImpl.cuh"
#include "LValueVectorImpl/LVectorBlock.cuh"
