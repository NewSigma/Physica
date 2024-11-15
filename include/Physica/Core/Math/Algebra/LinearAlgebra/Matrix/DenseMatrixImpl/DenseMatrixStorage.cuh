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

#include <Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.cuh>
#include "DenseMatrixStorage.h"

namespace Physica::Core {
    template<class T, size_t Row, size_t Col, class Allocator>
    class device_obj<DenseMatrixStorage<T, MatrixOption::Col | MatrixOption::Element, Row, Col, Allocator>>
            : public DenseMatrixDim<device_obj<DenseMatrixStorage<T, MatrixOption::Col | MatrixOption::Element, Row, Col, Allocator>>, Row, Col> {
        using host_obj = DenseMatrixStorage<T, MatrixOption::Col | MatrixOption::Element, Row, Col, Allocator>;
        using This = device_obj<host_obj>;
        using Dim = DenseMatrixDim<This, Row, Col>;
        using ArrayType = device_obj<typename std::conditional<is_scalar<T>::value, DenseVector<T, Row * Col>, Array<T, Row * Col>>::type>;
        using ElemType = typename ArrayType::ElemType;
    private:
        ArrayType arr;
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
        [[nodiscard]] host_obj toHost() const { return host_obj(arr.toHost(), {getRow(), getCol()}); }
        template<class... Args>
        __host__ __device__ void resize(size_t row, size_t col, Args&&... args);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        using Dim::getCol;
        using Dim::getRow;
        [[nodiscard]] __host__ __device__ ArrayType& asArray() noexcept { return arr; }
        [[nodiscard]] __host__ __device__ const ArrayType& asArray() const noexcept { return arr; }
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return arr.getLength(); }
        [[nodiscard]] __host__ __device__ ElemType* data_ptr(size_t row, size_t col) { return arr.data() + toIndex(row, col); }
        [[nodiscard]] __host__ __device__ const ElemType* data_ptr(size_t row, size_t col) const { return arr.data() + toIndex(row, col); }
    private:
        __host__ __device__ size_t toIndex(size_t r, size_t c) const { return getRow() * c + r; }
    };

    template<class T, size_t Row, size_t Col, class Allocator>
    class device_obj<DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Element, Row, Col, Allocator>>
            : public DenseMatrixDim<device_obj<DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Element, Row, Col, Allocator>>, Row, Col> {
        using host_obj = DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Element, Row, Col, Allocator>;
        using This = device_obj<host_obj>;
        using Dim = DenseMatrixDim<This, Row, Col>;
        using ArrayType = device_obj<typename std::conditional<is_scalar<T>::value, DenseVector<T, Row * Col>, Array<T, Row * Col>>::type>;
        using ElemType = typename ArrayType::ElemType;
    private:
        ArrayType arr;
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
        [[nodiscard]] host_obj toHost() const { return host_obj(arr.toHost(), {getRow(), getCol()}); }
        template<class... Args>
        __host__ __device__ void resize(size_t row, size_t col, Args&&... args);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        using Dim::getCol;
        using Dim::getRow;
        [[nodiscard]] __host__ __device__ ArrayType& asArray() noexcept { return arr; }
        [[nodiscard]] __host__ __device__ const ArrayType& asArray() const noexcept { return arr; }
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return arr.getLength(); }
        [[nodiscard]] __host__ __device__ ElemType* data_ptr(size_t row, size_t col) { return arr.data() + toIndex(row, col); }
        [[nodiscard]] __host__ __device__ const ElemType* data_ptr(size_t row, size_t col) const { return arr.data() + toIndex(row, col); }
    private:
        __host__ __device__ size_t toIndex(size_t r, size_t c) const { return getCol() * r + c; }
    };

    template<class T, size_t Row, size_t Col, class Allocator>
    class device_obj<DenseMatrixStorage<T, MatrixOption::Col | MatrixOption::Vector, Row, Col, Allocator>>
            : public DenseMatrixDim<device_obj<DenseMatrixStorage<T, MatrixOption::Col | MatrixOption::Vector, Row, Col, Allocator>>, Row, Col> {
        using host_obj = DenseMatrixStorage<T, MatrixOption::Col | MatrixOption::Vector, Row, Col, Allocator>;
        using This = device_obj<host_obj>;
        using Dim = DenseMatrixDim<This, Row, Col>;
        using VectorType = typename host_obj::VectorType;
        using ArrayType = device_obj<typename host_obj::ArrayType>;

        ArrayType array;
        size_t size;
    public:
        device_obj() = default;
        device_obj(size_t row, size_t col);
        device_obj(size_t row, size_t col, T value);
        device_obj(const host_obj& storage);
        device_obj(const This&) = default;
        device_obj(This&& obj) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] __device__ T& operator()(size_t r, size_t c);
        [[nodiscard]] __device__ const T& operator()(size_t r, size_t c) const;
        /* Operations */
        [[nodiscard]] host_obj toHost() const { return host_obj(array.toHost(), {getRow(), getCol()}); }
        void resize(size_t row, size_t col);
        inline void swap(This& __restrict obj) noexcept;
        /* Getters */
        using Dim::getCol;
        using Dim::getRow;
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return size; }
        [[nodiscard]] __host__ __device__ T* data_ptr(size_t row, size_t col);
        [[nodiscard]] __host__ __device__ const T* data_ptr(size_t row, size_t col) const;
    };

    template<class T, size_t Row, size_t Col, class Allocator>
    class device_obj<DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Vector, Row, Col, Allocator>>
            : public DenseMatrixDim<device_obj<DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Vector, Row, Col, Allocator>>, Row, Col> {
        using host_obj = DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Vector, Row, Col, Allocator>;
        using This = device_obj<host_obj>;
        using Dim = DenseMatrixDim<This, Row, Col>;
        using VectorType = typename host_obj::VectorType;
        using ArrayType = device_obj<typename host_obj::ArrayType>;

        ArrayType array;
        size_t size;
    public:
        device_obj() = default;
        device_obj(size_t row, size_t col);
        device_obj(size_t row, size_t col, T value);
        device_obj(const host_obj& storage);
        device_obj(const This&) = default;
        device_obj(This&& obj) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] __device__ T& operator()(size_t r, size_t c);
        [[nodiscard]] __device__ const T& operator()(size_t r, size_t c) const;
        /* Operations */
        [[nodiscard]] host_obj toHost() const { return host_obj(array.toHost(), {getRow(), getCol()}); }
        void resize(size_t row, size_t col);
        inline void swap(This& __restrict obj) noexcept;
        /* Getters */
        using Dim::getCol;
        using Dim::getRow;
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return size; }
        [[nodiscard]] __host__ __device__ T* data_ptr(size_t row, size_t col);
        [[nodiscard]] __host__ __device__ const T* data_ptr(size_t row, size_t col) const;
    };
}

#include "DenseMatrixStorageImpl.cuh"
