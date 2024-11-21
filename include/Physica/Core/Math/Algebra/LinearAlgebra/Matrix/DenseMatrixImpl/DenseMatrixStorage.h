/*
 * Copyright 2021-2024 Weibo He.
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

#include <cassert>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "DenseMatrixDim.h"

namespace Physica::Core {
    template<class T,
             int Option = MatrixOption::Col | MatrixOption::Vector,
             size_t Row = Dynamic,
             size_t Col = Dynamic,
             class Allocator = HostAllocator<T>>
    class DenseMatrixStorage;

    template<class T, size_t Row, size_t Col, class Allocator>
    class DenseMatrixStorage<T, MatrixOption::Col | MatrixOption::Element, Row, Col, Allocator>
            : public DenseMatrixDim<DenseMatrixStorage<T, MatrixOption::Col | MatrixOption::Element, Row, Col, Allocator>, Row, Col> {
        using This = DenseMatrixStorage<T, MatrixOption::Col | MatrixOption::Element, Row, Col, Allocator>;
        using Dim = DenseMatrixDim<This, Row, Col>;
        using RefTy = T::RefTy;
        using ConstRefTy = T::ConstRefTy;
        using ArrayType = std::conditional<Scalar<T>, DenseVector<T, Row * Col>, Array<T, Row * Col>>::type;
    public:
        using InitializerType = T;
    private:
        ArrayType arr;
    public:
        DenseMatrixStorage() = default;
        DenseMatrixStorage(size_t row, size_t col) : Dim(row, col), arr(row * col) {}
        DenseMatrixStorage(size_t row, size_t col, T value) : Dim(row, col), arr(row * col, std::move(value)) {}
        DenseMatrixStorage(std::initializer_list<T> list) : arr(list) {
            assert(Row != Dynamic && Col != Dynamic && "[Error]: Row or Col is unknown");
        }
        DenseMatrixStorage(const This&) = default;
        DenseMatrixStorage(This&&) noexcept = default;
        ~DenseMatrixStorage() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] RefTy operator()(size_t r, size_t c) {
            assert(r < getRow() && c < getCol());
            return arr[toIndex(r, c)];
        }

        [[nodiscard]] ConstRefTy operator()(size_t r, size_t c) const {
            assert(r < getRow() && c < getCol());
            return arr[toIndex(r, c)];
        }
        /* Operations */
        template<class... Args>
        void resize(size_t row, size_t col, Args&&... args);
        [[nodiscard]] inline auto toDevice() const;
        [[nodiscard]] inline auto toDeviceAsync() const;
        inline void toDevice(device_obj<This>& obj) const;
        inline void toDeviceAsync(device_obj<This>& obj) const;

        void swap(This& __restrict obj) noexcept;
        void swap_row(size_t r1, size_t r2) noexcept;
        void swap_col(size_t c1, size_t r1) noexcept;
        /* Getters */
        using Dim::getCol;
        using Dim::getRow;
        [[nodiscard]] ArrayType& asArray() noexcept { return arr; }
        [[nodiscard]] const ArrayType& asArray() const noexcept { return arr; }
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return arr.getLength(); }
        [[nodiscard]] __host__ __device__ T* data_ptr(size_t row, size_t col) { return arr.data() + toIndex(row, col); }
        [[nodiscard]] __host__ __device__ const T* data_ptr(size_t row, size_t col) const { return arr.data() + toIndex(row, col); }
    private:
        DenseMatrixStorage(ArrayType arr_, Dim dim) : Dim(std::move(dim)), arr(std::move(arr_)) {}
        __host__ __device__ size_t toIndex(size_t r, size_t c) const { return getRow() * c + r; }
        /* Friends */
        friend class device_obj<This>;
    };

    template<class T, size_t Row, size_t Col, class Allocator>
    class DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Element, Row, Col, Allocator>
            : public DenseMatrixDim<DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Element, Row, Col, Allocator>, Row, Col> {
        using This = DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Element, Row, Col, Allocator>;
        using Dim = DenseMatrixDim<This, Row, Col>;
        using RefTy = T::RefTy;
        using ConstRefTy = T::ConstRefTy;
        using ArrayType = std::conditional<Scalar<T>, DenseVector<T, Row * Col>, Array<T, Row * Col>>::type;
    public:
        using InitializerType = T;
    private:
        ArrayType arr;
    public:
        DenseMatrixStorage() = default;
        DenseMatrixStorage(size_t row, size_t col) : Dim(row, col), arr(row * col) {}
        DenseMatrixStorage(size_t row, size_t col, T value) : Dim(row, col), arr(row * col, std::move(value)) {}
        DenseMatrixStorage(std::initializer_list<T> list) : arr(list) {}
        DenseMatrixStorage(const This&) = default;
        DenseMatrixStorage(This&&) noexcept = default;
        ~DenseMatrixStorage() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] RefTy& operator()(size_t r, size_t c) {
            assert(r < getRow() && c < getCol());
            return arr[toIndex(r, c)];
        }

        [[nodiscard]] ConstRefTy operator()(size_t r, size_t c) const {
            assert(r < getRow() && c < getCol());
            return arr[toIndex(r, c)];
        }
        /* Operations */
        template<class... Args>
        void resize(size_t row, size_t col, Args&&... args);
        [[nodiscard]] inline auto toDevice() const;
        [[nodiscard]] inline auto toDeviceAsync() const;
        inline void toDevice(device_obj<This>& obj) const;
        inline void toDeviceAsync(device_obj<This>& obj) const;

        void swap(This& __restrict obj) noexcept;
        void swap_row(size_t r1, size_t r2) noexcept;
        void swap_col(size_t c1, size_t r1) noexcept;
        /* Getters */
        using Dim::getCol;
        using Dim::getRow;
        [[nodiscard]] ArrayType& asArray() noexcept { return arr; }
        [[nodiscard]] const ArrayType& asArray() const noexcept { return arr; }
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return arr.getLength(); }
        [[nodiscard]] __host__ __device__ T* data_ptr(size_t row, size_t col) { return arr.data() + toIndex(row, col); }
        [[nodiscard]] __host__ __device__ const T* data_ptr(size_t row, size_t col) const { return arr.data() + toIndex(row, col); }
    private:
        DenseMatrixStorage(ArrayType arr_, Dim dim) : Dim(std::move(dim)), arr(std::move(arr_)) {}
        __host__ __device__ size_t toIndex(size_t r, size_t c) const { return getCol() * r + c; }
        /* Friends */
        friend class device_obj<This>;
    };

    template<class T, size_t Row, size_t Col, class Allocator>
    class DenseMatrixStorage<T, MatrixOption::Col | MatrixOption::Vector, Row, Col, Allocator>
            : public DenseMatrixDim<DenseMatrixStorage<T, MatrixOption::Col | MatrixOption::Vector, Row, Col, Allocator>, Row, Col> {
        using This = DenseMatrixStorage<T, MatrixOption::Col | MatrixOption::Vector, Row, Col, Allocator>;
        using Dim = DenseMatrixDim<This, Row, Col>;
        using RefTy = T::RefTy;
        using ConstRefTy = T::ConstRefTy;
    public:
        using VectorType = std::conditional<Scalar<T>, DenseVector<T, Row>, Array<T, Row>>::type;
        using InitializerType = VectorType;
    private:
        using AllocatorV = Allocator::template rebind_alloc<VectorType>;
        using ArrayType = Array<VectorType, Col, AllocatorV>;

        ArrayType array;
    public:
        DenseMatrixStorage() = default;
        DenseMatrixStorage(size_t row, size_t col) : Dim(row, col), array(col, row) {}
        DenseMatrixStorage(size_t row, size_t col, T value) : Dim(row, col), array(col, row, std::move(value)) {}
        DenseMatrixStorage(std::initializer_list<VectorType> list) : array(list) {
            const size_t length = array.getLength();
            const size_t size = getSize();
            Dim::resize(size / length, length);
        }
        DenseMatrixStorage(const This&) = default;
        DenseMatrixStorage(This&&) noexcept = default;
        ~DenseMatrixStorage() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] RefTy operator()(size_t r, size_t c) {
            assert(r < getRow() && c < getCol());
            return array[c][r];
        }

        [[nodiscard]] ConstRefTy operator()(size_t r, size_t c) const {
            assert(r < getRow() && c < getCol());
            return array[c][r];
        }
        /* Operations */
        template<class... Args>
        void resize(size_t row, size_t col, Args&&... args);
        [[nodiscard]] inline auto toDevice() const;

        void swap(This& __restrict obj) noexcept;
        void swap_row(size_t r1, size_t r2) noexcept;
        void swap_col(size_t c1, size_t r1) noexcept;
        /* Getters */
        using Dim::getCol;
        using Dim::getRow;
        [[nodiscard]] ArrayType& asArray() noexcept { return array; }
        [[nodiscard]] const ArrayType& asArray() const noexcept { return array; }
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept {
            return array.getLength() == 0 ? 0 : array.getLength() * array[0].getLength();
        }
        [[nodiscard]] __host__ __device__ T* data_ptr(size_t row, size_t col) { return array[col].data() + row; }
        [[nodiscard]] __host__ __device__ const T* data_ptr(size_t row, size_t col) const { return array[col].data() + row; }
    private:
        DenseMatrixStorage(ArrayType array, Dim dim) : Dim(std::move(dim)), array(std::move(array)) {}
        friend class device_obj<This>;
    };

    template<class T, size_t Row, size_t Col, class Allocator>
    class DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Vector, Row, Col, Allocator>
            : public DenseMatrixDim<DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Vector, Row, Col, Allocator>, Row, Col> {
        using This = DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Vector, Row, Col, Allocator>;
        using Dim = DenseMatrixDim<This, Row, Col>;
        using RefTy = T::RefTy;
        using ConstRefTy = T::ConstRefTy;
    public:
        using VectorType = std::conditional<Scalar<T>, DenseVector<T, Col>, Array<T, Col>>::type;
        using InitializerType = VectorType;
    private:
        using AllocatorV = Allocator::template rebind_alloc<VectorType>;
        using ArrayType = Array<VectorType, Row, AllocatorV>;

        ArrayType array;
    public:
        DenseMatrixStorage() = default;
        DenseMatrixStorage(size_t row, size_t col) : Dim(row, col), array(row, col) {}
        DenseMatrixStorage(size_t row, size_t col, T value) : Dim(row, col), array(row, col, std::move(value)) {}
        DenseMatrixStorage(std::initializer_list<VectorType> list) : array(list) {
            const size_t length = array.getLength();
            const size_t size = getSize();
            Dim::resize(length, size / length);
        }
        DenseMatrixStorage(const This&) = default;
        DenseMatrixStorage(This&&) noexcept = default;
        ~DenseMatrixStorage() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] RefTy operator()(size_t r, size_t c) {
            assert(r < getRow() && c < getCol());
            return array[r][c];
        }

        [[nodiscard]] ConstRefTy operator()(size_t r, size_t c) const {
            assert(r < getRow() && c < getCol());
            return array[r][c];
        }
        /* Operations */
        template<class... Args>
        void resize(size_t row, size_t col, Args&&... args);
        [[nodiscard]] inline auto toDevice() const;

        void swap(This& __restrict obj) noexcept;
        void swap_row(size_t r1, size_t r2) noexcept;
        void swap_col(size_t c1, size_t r1) noexcept;
        /* Getters */
        using Dim::getCol;
        using Dim::getRow;
        [[nodiscard]] ArrayType& asArray() noexcept { return array; }
        [[nodiscard]] const ArrayType& asArray() const noexcept { return array; }
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept {
            return array.getLength() == 0 ? 0 : array.getLength() * array[0].getLength();
        }
        [[nodiscard]] __host__ __device__ T* data_ptr(size_t row, size_t col) { return array[row].data() + col; }
        [[nodiscard]] __host__ __device__ const T* data_ptr(size_t row, size_t col) const { return array[row].data() + col; }
    private:
        DenseMatrixStorage(ArrayType array, Dim dim) : Dim(std::move(dim)), array(std::move(array)) {}
        friend class device_obj<This>;
    };
}

#include "DenseMatrixStorageImpl.h"
