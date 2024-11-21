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

#include "RValueMatrix.cuh"

namespace Physica::Core {
    template<class Derived>
    class device_obj<LValueMatrix<Derived>> : public device_obj<RValueMatrix<Derived>> {
        using Base = device_obj<RValueMatrix<Derived>>;
    public:
        using typename Base::ScalarType;
        using Base::RowAtCompile;
        using Base::ColAtCompile;
    protected:
        using PtrTy = ScalarType::PtrTy;
        using ConstPtrTy = ScalarType::ConstPtrTy;
    public:
        /* Operators */
        device_obj& operator=(const device_obj& m) = delete;
        device_obj& operator=(device_obj&& m) = delete;
        template<Matrix M>
        __host__ __device__ device_obj<Derived>& operator=(const device_obj<M>& m);
        __device__ device_obj<Derived>& operator=(const ScalarType& x);
        [[nodiscard]] __device__ ScalarType& operator()(size_t row, size_t col) { return *data_ptr(row, col); }
        [[nodiscard]] __device__ const ScalarType& operator()(size_t row, size_t col) const { return *data_ptr(row, col); }
        __device__ void operator+=(const ScalarType& s) { Base::getDerived() = Base::getDerived() + s; }
        __device__ void operator-=(const ScalarType& s) { Base::getDerived() = Base::getDerived() - s; }
        __device__ void operator*=(const ScalarType& s) { Base::getDerived() = Base::getDerived() * s; }
        __device__ void operator/=(const ScalarType& s) { Base::getDerived() = Base::getDerived() / s; }
        template<Matrix M> __device__ void operator+=(const device_obj<M>& m) { Base::getDerived() = Base::getDerived() + m; }
        template<Matrix M> __device__ void operator-=(const device_obj<M>& m) { Base::getDerived() = Base::getDerived() - m; }
        /* Getters */
        template<Side Owner = GetSide()>
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const { return *data_ptr(row, col); }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr(size_t row, size_t col) { return Base::getDerived().data_ptr(row, col); }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr(size_t row, size_t col) const { return Base::getDerived().data_ptr(row, col); }
        [[nodiscard]] __device__ inline ScalarType& refFromMajorMinor(size_t major, size_t minor);
        [[nodiscard]] __device__ inline const ScalarType& refFromMajorMinor(size_t major, size_t minor) const;
        [[nodiscard]] __host__ __device__ device_obj<LValueFlatten<Derived>> flatten() { return {*this}; }
        [[nodiscard]] __host__ __device__ const device_obj<LValueFlatten<Derived>> flatten() const { return {*this}; }
    protected:
        device_obj() = default;
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
    };
}

#include "LValueMatrixImpl/LValueMatrixImpl.cuh"
#include "LValueMatrixImpl/LValueFlatten.cuh"
