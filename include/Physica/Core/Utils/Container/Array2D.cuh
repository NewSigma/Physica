/*
 * Copyright 2022-2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.cuh" // IWYU pragma: export
#include "Array2D.h"

namespace Physica {
    template<class T, int Option, size_t Row, size_t Col, class Allocator>
    class device_obj<Array2D<T, Option, Row, Col, Allocator>> {
        using host_obj = Array2D<T, Option, Row, Col, Allocator>;
        using This = device_obj<host_obj>;
        using ArrayType = device_obj<typename host_obj::ArrayType>;
        using IndexType = host_obj::IndexType;

        constexpr static bool isVectorStorage = host_obj::isVectorStorage;
        constexpr static bool isDynamicArray = host_obj::isDynamicArray;
        constexpr static bool isColMajor = host_obj::isColMajor;
        using ElemType = std::conditional<isVectorStorage, T, typename ArrayType::ElemType>::type;
        using SizeType = std::conditional<isVectorStorage, size_t, PlainStruct<void>>::type;
    private:
        ArrayType arr;
        [[no_unique_address]] IndexType r = 0;
        SizeType size = 0; // Double [[no_unique_address]] makes NVCC 12.8, arch 75 miscompilation
    public:
        device_obj() = default;
        __host__ __device__ device_obj(size_t row, size_t col);
        __host__ __device__ device_obj(size_t row, size_t col, ElemType value);
        device_obj(const host_obj& storage);
        device_obj(const This&) = default;
        device_obj(This&& obj) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] __device__ ElemType& operator()(size_t r, size_t c);
        [[nodiscard]] __device__ const ElemType& operator()(size_t r, size_t c) const;
        /* Operations */
        template<class... Args>
        __host__ __device__ void resize(size_t row, size_t col, Args&&... args);

        [[nodiscard]] host_obj toHost() const;
        [[nodiscard]] host_obj toHostAsync() const;
        void toHost(host_obj& obj) const;
        void toHostAsync(host_obj& obj) const;

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept;
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept;
        [[nodiscard]] __host__ __device__ ArrayType& asArray() noexcept { return arr; }
        [[nodiscard]] __host__ __device__ const ArrayType& asArray() const noexcept { return arr; }
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept;
        [[nodiscard]] __host__ __device__ ElemType* data_ptr(size_t row, size_t col);
        [[nodiscard]] __host__ __device__ const ElemType* data_ptr(size_t row, size_t col) const;
    private:
        __host__ __device__ size_t toIndex1D(size_t r, size_t c) const;
    };
}

#include "ArrayImpl/Array2DImpl.cuh"
