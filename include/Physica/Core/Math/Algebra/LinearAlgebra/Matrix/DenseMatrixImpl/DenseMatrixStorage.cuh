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

#include <Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.cuh>
#include "DenseMatrixStorage.h"

namespace Physica::Core {
    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    class device_obj<DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>>
            : public Utils::device_obj<Utils::Array<T, Row * Column, MaxRow * MaxColumn>>
            , public DenseMatrixDim<device_obj<DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>>, Row, Column, MaxRow, MaxColumn> {
        using host_obj = DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>;
        using This = device_obj<host_obj>;
        using Base = Utils::device_obj<Utils::Array<T, Row * Column, MaxRow * MaxColumn>>;
        using Dim = DenseMatrixDim<This, Row, Column, MaxRow, MaxColumn>;
    public:
        device_obj() = default;
        __host__ __device__ device_obj(size_t row, size_t column);
        __host__ __device__ device_obj(size_t row, size_t column, T value);
        device_obj(const host_obj& storage);
        device_obj(const This&) = default;
        device_obj(This&& obj) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] __device__ T& operator()(size_t r, size_t c);
        [[nodiscard]] __device__ const T& operator()(size_t r, size_t c) const;
        /* Operations */
        [[nodiscard]] host_obj toHost() const { return host_obj(Base::toHost(), {getRow(), getColumn()}); }
        template<class... Args>
        __host__ __device__ void resize(size_t row, size_t column, Args&&... args);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        using Dim::getColumn;
        using Dim::getRow;
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return Base::getLength(); }
        [[nodiscard]] __host__ __device__ T* data_ptr(size_t row, size_t column) { return Base::data() + toIndex(row, column); }
        [[nodiscard]] __host__ __device__ const T* data_ptr(size_t row, size_t column) const { return Base::data() + toIndex(row, column); }
    private:
        __host__ __device__ size_t toIndex(size_t r, size_t c) const { return getRow() * c + r; }

        using Base::getLength;
    };

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    class device_obj<DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>>
            : public Utils::device_obj<Utils::Array<T, Row * Column, MaxRow * MaxColumn>>
            , public DenseMatrixDim<device_obj<DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>>, Row, Column, MaxRow, MaxColumn> {
        using host_obj = DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>;
        using This = device_obj<host_obj>;
        using Base = Utils::device_obj<Utils::Array<T, Row * Column, MaxRow * MaxColumn>>;
        using Dim = DenseMatrixDim<This, Row, Column, MaxRow, MaxColumn>;
    public:
        device_obj() = default;
        __host__ __device__ device_obj(size_t row, size_t column);
        __host__ __device__ device_obj(size_t row, size_t column, T value);
        device_obj(const host_obj& storage);
        device_obj(const This&) = default;
        device_obj(This&& obj) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] __device__ T& operator()(size_t r, size_t c);
        [[nodiscard]] __device__ const T& operator()(size_t r, size_t c) const;
        /* Operations */
        [[nodiscard]] host_obj toHost() const { return host_obj(Base::toHost(), {getRow(), getColumn()}); }
        template<class... Args>
        __host__ __device__ void resize(size_t row, size_t column, Args&&... args);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        using Dim::getColumn;
        using Dim::getRow;
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return Base::getLength(); }
        [[nodiscard]] __host__ __device__ T* data_ptr(size_t row, size_t column) { return Base::data() + toIndex(row, column); }
        [[nodiscard]] __host__ __device__ const T* data_ptr(size_t row, size_t column) const { return Base::data() + toIndex(row, column); }
    private:
        __host__ __device__ size_t toIndex(size_t r, size_t c) const { return getColumn() * r + c; }

        using Base::getLength;
    };

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    class device_obj<DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>>
            : public DenseMatrixDim<device_obj<DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>>, Row, Column, MaxRow, MaxColumn> {
        using host_obj = DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>;
        using This = device_obj<host_obj>;
        using Dim = DenseMatrixDim<This, Row, Column, MaxRow, MaxColumn>;
        using VectorType = typename host_obj::VectorType;
        using ArrayType = Physica::Utils::device_obj<typename host_obj::ArrayType>;

        ArrayType array;
        size_t size;
    public:
        device_obj() = default;
        device_obj(size_t row, size_t column);
        device_obj(size_t row, size_t column, T value);
        device_obj(const host_obj& storage);
        device_obj(const This&) = default;
        device_obj(This&& obj) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] __device__ T& operator()(size_t r, size_t c);
        [[nodiscard]] __device__ const T& operator()(size_t r, size_t c) const;
        /* Operations */
        [[nodiscard]] host_obj toHost() const { return host_obj(array.toHost(), {getRow(), getColumn()}); }
        void resize(size_t row, size_t column);
        inline void swap(This& __restrict obj) noexcept;
        /* Getters */
        using Dim::getColumn;
        using Dim::getRow;
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return size; }
        [[nodiscard]] __host__ __device__ T* data_ptr(size_t row, size_t column);
        [[nodiscard]] __host__ __device__ const T* data_ptr(size_t row, size_t column) const;
    };

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    class device_obj<DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>>
            : public DenseMatrixDim<device_obj<DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>>, Row, Column, MaxRow, MaxColumn> {
        using host_obj = DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>;
        using This = device_obj<host_obj>;
        using Dim = DenseMatrixDim<This, Row, Column, MaxRow, MaxColumn>;
        using VectorType = typename host_obj::VectorType;
        using ArrayType = Physica::Utils::device_obj<typename host_obj::ArrayType>;

        ArrayType array;
        size_t size;
    public:
        device_obj() = default;
        device_obj(size_t row, size_t column);
        device_obj(size_t row, size_t column, T value);
        device_obj(const host_obj& storage);
        device_obj(const This&) = default;
        device_obj(This&& obj) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] __device__ T& operator()(size_t r, size_t c);
        [[nodiscard]] __device__ const T& operator()(size_t r, size_t c) const;
        /* Operations */
        [[nodiscard]] host_obj toHost() const { return host_obj(array.toHost(), {getRow(), getColumn()}); }
        void resize(size_t row, size_t column);
        inline void swap(This& __restrict obj) noexcept;
        /* Getters */
        using Dim::getColumn;
        using Dim::getRow;
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return size; }
        [[nodiscard]] __host__ __device__ T* data_ptr(size_t row, size_t column);
        [[nodiscard]] __host__ __device__ const T* data_ptr(size_t row, size_t column) const;
    };
}

#include "DenseMatrixStorageImpl.cuh"
