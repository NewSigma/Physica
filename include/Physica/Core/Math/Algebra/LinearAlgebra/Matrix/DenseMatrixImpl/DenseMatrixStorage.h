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
#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixOption.h>
#include <Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h>
#include "DenseMatrixDim.h"

namespace Physica::Core {
    template<class T,
             int Option = MatrixOption::Column | MatrixOption::Vector,
             size_t Row = Dynamic,
             size_t Column = Dynamic,
             class Allocator = HostAllocator<T>>
    class DenseMatrixStorage;

    template<class T, size_t Row, size_t Column, class Allocator>
    class DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Element, Row, Column, Allocator>
            : public Array<T, Row * Column, Allocator>
            , public DenseMatrixDim<DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Element, Row, Column, Allocator>, Row, Column> {
        using This = DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Element, Row, Column, Allocator>;
        using Dim = DenseMatrixDim<This, Row, Column>;
    public:
        using Base = Array<T, Row * Column, Allocator>;
        using InitializerType = T;
    public:
        DenseMatrixStorage() = default;
        DenseMatrixStorage(size_t row, size_t column) : Base(row * column), Dim(row, column) {}
        DenseMatrixStorage(size_t row, size_t column, T value) : Base(row * column, std::move(value)), Dim(row, column) {}
        DenseMatrixStorage(std::initializer_list<T> list) : Base(list) {
            assert(Row != Dynamic && Column != Dynamic && "[Error]: Row or Column is unknown");
        }
        DenseMatrixStorage(const This&) = default;
        DenseMatrixStorage(This&&) noexcept = default;
        ~DenseMatrixStorage() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] T& operator()(size_t r, size_t c) {
            assert(r < getRow() && c < getColumn());
            return Base::operator[](toIndex(r, c));
        }

        [[nodiscard]] const T& operator()(size_t r, size_t c) const {
            assert(r < getRow() && c < getColumn());
            return Base::operator[](toIndex(r, c));
        }
        /* Operations */
        template<class... Args>
        void resize(size_t row, size_t column, Args&&... args);
        void rowSwap(size_t r1, size_t r2);
        void columnSwap(size_t c1, size_t r1);
        [[nodiscard]] inline auto toDevice() const;
        [[nodiscard]] inline auto toDeviceAsync() const;
        inline void toDevice(device_obj<This>& obj) const;
        inline void toDeviceAsync(device_obj<This>& obj) const;
        void swap(This& __restrict obj) noexcept {
            assert(this != &obj && "[Error]: Self swap is likely a bug");
            Base::swap(obj);
            Dim::swap(obj);
        }
        /* Getters */
        using Dim::getColumn;
        using Dim::getRow;
        [[nodiscard]] Base& asArray() noexcept { return *this; }
        [[nodiscard]] const Base& asArray() const noexcept { return *this; }
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return Base::getLength(); }
        [[nodiscard]] __host__ __device__ T* data_ptr(size_t row, size_t column) { return Base::data() + toIndex(row, column); }
        [[nodiscard]] __host__ __device__ const T* data_ptr(size_t row, size_t column) const { return Base::data() + toIndex(row, column); }
    private:
        DenseMatrixStorage(Base array, Dim dim) : Base(std::move(array)), Dim(std::move(dim)) {}
        __host__ __device__ size_t toIndex(size_t r, size_t c) const { return getRow() * c + r; }

        using Base::getCapacity;
        using Base::getLength;
        friend class device_obj<This>;
    };

    template<class T, size_t Row, size_t Column, class Allocator>
    class DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Element, Row, Column, Allocator>
            : public Array<T, Row * Column, Allocator>
            , public DenseMatrixDim<DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Element, Row, Column, Allocator>, Row, Column> {
        using This = DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Element, Row, Column, Allocator>;
        using Dim = DenseMatrixDim<This, Row, Column>;
    public:
        using Base = Array<T, Row * Column, Allocator>;
        using InitializerType = T;
    public:
        DenseMatrixStorage() = default;
        DenseMatrixStorage(size_t row, size_t column) : Base(row * column), Dim(row, column) {}
        DenseMatrixStorage(size_t row, size_t column, T value) : Base(row * column, std::move(value)), Dim(row, column) {}
        DenseMatrixStorage(std::initializer_list<T> list) : Base(list) {}
        DenseMatrixStorage(const This&) = default;
        DenseMatrixStorage(This&&) noexcept = default;
        ~DenseMatrixStorage() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] T& operator()(size_t r, size_t c) {
            assert(r < getRow() && c < getColumn());
            return Base::operator[](toIndex(r, c));
        }

        [[nodiscard]] const T& operator()(size_t r, size_t c) const {
            assert(r < getRow() && c < getColumn());
            return Base::operator[](toIndex(r, c));
        }
        /* Operations */
        template<class... Args>
        void resize(size_t row, size_t column, Args&&... args);
        void rowSwap(size_t r1, size_t r2);
        void columnSwap(size_t c1, size_t r1);
        [[nodiscard]] inline auto toDevice() const;
        [[nodiscard]] inline auto toDeviceAsync() const;
        inline void toDevice(device_obj<This>& obj) const;
        inline void toDeviceAsync(device_obj<This>& obj) const;
        void swap(This& __restrict obj) noexcept {
            assert(this != &obj && "[Error]: Self swap is likely a bug");
            Base::swap(obj);
            Dim::swap(obj);
        }
        /* Getters */
        using Dim::getColumn;
        using Dim::getRow;
        [[nodiscard]] Base& asArray() noexcept { return *this; }
        [[nodiscard]] const Base& asArray() const noexcept { return *this; }
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return Base::getLength(); }
        [[nodiscard]] __host__ __device__ T* data_ptr(size_t row, size_t column) { return Base::data() + toIndex(row, column); }
        [[nodiscard]] __host__ __device__ const T* data_ptr(size_t row, size_t column) const { return Base::data() + toIndex(row, column); }
    private:
        DenseMatrixStorage(Base array, Dim dim) : Base(std::move(array)), Dim(std::move(dim)) {}
        __host__ __device__ size_t toIndex(size_t r, size_t c) const { return getColumn() * r + c; }

        using Base::getCapacity;
        using Base::getLength;
        friend class device_obj<This>;
    };

    template<class T, size_t Row, size_t Column, class Allocator>
    class DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Vector, Row, Column, Allocator>
            : public DenseMatrixDim<DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Vector, Row, Column, Allocator>, Row, Column> {
        using This = DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Vector, Row, Column, Allocator>;
        using Dim = DenseMatrixDim<This, Row, Column>;
    public:
        using VectorType = typename std::conditional<is_scalar<T>::value, Vector<T, Row, Allocator>, Array<T, Row, Allocator>>::type;
        using InitializerType = VectorType;
    private:
        using AllocatorV = typename ChangeAllocatorValueType<Allocator, VectorType>::Type;
        using ArrayType = Array<VectorType, Column, AllocatorV>;

        ArrayType array;
    public:
        DenseMatrixStorage() = default;
        DenseMatrixStorage(size_t row, size_t column) : Dim(row, column), array(column, row) {}
        DenseMatrixStorage(size_t row, size_t column, T value) : Dim(row, column), array(column, row, std::move(value)) {}
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
        [[nodiscard]] T& operator()(size_t r, size_t c) {
            assert(r < getRow() && c < getColumn());
            return array[c][r];
        }

        [[nodiscard]] const T& operator()(size_t r, size_t c) const {
            assert(r < getRow() && c < getColumn());
            return array[c][r];
        }
        /* Operations */
        template<class... Args>
        void resize(size_t row, size_t column, Args&&... args);
        void rowSwap(size_t r1, size_t r2);
        void columnSwap(size_t c1, size_t r1);
        [[nodiscard]] inline auto toDevice() const;
        void swap(This& __restrict obj) noexcept {
            assert(this != &obj && "[Error]: Self swap is likely a bug");
            array.swap(obj.array);
            Dim::swap(obj);
        }
        /* Getters */
        using Dim::getColumn;
        using Dim::getRow;
        [[nodiscard]] ArrayType& asArray() noexcept { return array; }
        [[nodiscard]] const ArrayType& asArray() const noexcept { return array; }
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept {
            return array.getLength() == 0 ? 0 : array.getLength() * array[0].getLength();
        }
        [[nodiscard]] __host__ __device__ T* data_ptr(size_t row, size_t column) { return array[column].data() + row; }
        [[nodiscard]] __host__ __device__ const T* data_ptr(size_t row, size_t column) const { return array[column].data() + row; }
    private:
        DenseMatrixStorage(ArrayType array, Dim dim) : Dim(std::move(dim)), array(std::move(array)) {}
        friend class device_obj<This>;
    };

    template<class T, size_t Row, size_t Column, class Allocator>
    class DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Vector, Row, Column, Allocator>
            : public DenseMatrixDim<DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Vector, Row, Column, Allocator>, Row, Column> {
        using This = DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Vector, Row, Column, Allocator>;
        using Dim = DenseMatrixDim<This, Row, Column>;
    public:
        using VectorType = typename std::conditional<is_scalar<T>::value, Vector<T, Column, Allocator>, Array<T, Column, Allocator>>::type;
        using InitializerType = VectorType;
    private:
        using AllocatorV = typename ChangeAllocatorValueType<Allocator, VectorType>::Type;
        using ArrayType = Array<VectorType, Row, AllocatorV>;

        ArrayType array;
    public:
        DenseMatrixStorage() = default;
        DenseMatrixStorage(size_t row, size_t column) : Dim(row, column), array(row, column) {}
        DenseMatrixStorage(size_t row, size_t column, T value) : Dim(row, column), array(row, column, std::move(value)) {}
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
        [[nodiscard]] T& operator()(size_t r, size_t c) {
            assert(r < getRow() && c < getColumn());
            return array[r][c];
        }

        [[nodiscard]] const T& operator()(size_t r, size_t c) const {
            assert(r < getRow() && c < getColumn());
            return array[r][c];
        }
        /* Operations */
        template<class... Args>
        void resize(size_t row, size_t column, Args&&... args);
        void rowSwap(size_t r1, size_t r2);
        void columnSwap(size_t c1, size_t r1);
        [[nodiscard]] inline auto toDevice() const;
        void swap(This& __restrict obj) noexcept {
            assert(this != &obj && "[Error]: Self swap is likely a bug");
            array.swap(obj.array);
            Dim::swap(obj);
        }
        /* Getters */
        using Dim::getColumn;
        using Dim::getRow;
        [[nodiscard]] ArrayType& asArray() noexcept { return array; }
        [[nodiscard]] const ArrayType& asArray() const noexcept { return array; }
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept {
            return array.getLength() == 0 ? 0 : array.getLength() * array[0].getLength();
        }
        [[nodiscard]] __host__ __device__ T* data_ptr(size_t row, size_t column) { return array[row].data() + column; }
        [[nodiscard]] __host__ __device__ const T* data_ptr(size_t row, size_t column) const { return array[row].data() + column; }
    private:
        DenseMatrixStorage(ArrayType array, Dim dim) : Dim(std::move(dim)), array(std::move(array)) {}
        friend class device_obj<This>;
    };
}

#include "DenseMatrixStorageImpl.h"
