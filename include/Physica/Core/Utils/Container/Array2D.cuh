/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.cuh" // IWYU pragma: export
#include "Array2D.h"

namespace Physica {
    template<class T, int Major, size_t Row, size_t Col, class Allocator>
    class device_obj<Array2D<T, Major, Row, Col, Allocator>> {
        using host_obj = Array2D<T, Major, Row, Col, Allocator>;
        using This = device_obj<host_obj>;
        using ArrayType = device_obj<typename host_obj::ArrayType>;
        using IndexType = host_obj::IndexType;

        constexpr static bool isDynamicArray = host_obj::isDynamicArray;
        constexpr static bool isColMajor = host_obj::isColMajor;
    private:
        ArrayType arr;
        [[no_unique_address]] IndexType r = 0;
    public:
        device_obj() = default;
        explicit __host__ __device__ device_obj(size_t order);
        __host__ __device__ device_obj(size_t row, size_t col);
        __host__ __device__ device_obj(size_t row, size_t col, T value);
        device_obj(const host_obj& storage);
        device_obj(const This&) = default;
        device_obj(This&& obj) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] __device__ T& operator[](size_t r, size_t c);
        [[nodiscard]] __device__ const T& operator[](size_t r, size_t c) const;
        /* Operations */
        __host__ __device__ void resize(size_t row, size_t col, auto&&... args);
        __host__ __device__ void resize(size_t order);
        void reserve(size_t size);

        [[nodiscard]] host_obj toHost() const;
        [[nodiscard]] host_obj toHostAsync() const;
        void toHost(host_obj& obj) const;
        void toHostAsync(host_obj& obj) const;

        void zeros();
        void junk();
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ auto&& asArray(this auto&& self) noexcept { return self.arr; }
        [[nodiscard, gnu::returns_nonnull]] __host__ __device__ auto* data(this auto&&) noexcept;
        [[nodiscard, gnu::returns_nonnull]] __host__ __device__ auto* data_ptr(this auto&&, size_t row, size_t col) noexcept;
        [[nodiscard]] __host__ __device__ constexpr size_t getRow() const noexcept;
        [[nodiscard]] __host__ __device__ constexpr size_t getCol() const noexcept;
        [[nodiscard]] __host__ __device__ constexpr size_t getSize() const noexcept;
        [[nodiscard]] __host__ __device__ constexpr size_t getOrder() const noexcept;
    private:
        [[nodiscard]] __host__ __device__ size_t toIndex1D(size_t r, size_t c) const;
    };
}

#include "ArrayImpl/Array2DImpl.cuh"
