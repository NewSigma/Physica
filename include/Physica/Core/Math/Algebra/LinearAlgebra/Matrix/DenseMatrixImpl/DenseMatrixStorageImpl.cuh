/*
 * Copyright 2024 Weibo He.
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

namespace Physica::Core {
    //////////////////////////////////////////////Column-Element//////////////////////////////////////////////
    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    __host__ __device__ device_obj<DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>>::
            device_obj(size_t row, size_t column) {
        resize(row, column);
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    __host__ __device__ device_obj<DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>>::
            device_obj(size_t row, size_t column, T value) {
        resize(row, column, std::move(value));
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    device_obj<DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>>::
            device_obj(const host_obj& storage) : Base(storage), Dim(storage.getRow(), storage.getColumn()) {}

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    __device__ T& device_obj<DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>>::
            operator()(size_t r, size_t c) {
        assert(r < getRow() && c < getColumn());
        return Base::operator[](toIndex(r, c));
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    __device__ const T&
    device_obj<DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>>::
            operator()(size_t r, size_t c) const {
        assert(r < getRow() && c < getColumn());
        return Base::operator[](toIndex(r, c));
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    template<class... Args>
    __host__ __device__ void device_obj<DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>>::
            resize(size_t row, size_t column, Args&&... args) {
    #ifdef __CUDA_ARCH__
        assert(Row * Column != Dynamic && "[Error]: Do not allocate dynamic matrix in device code");
    #endif
        Base::resize(row * column, std::forward<Args>(args)...);
        Dim::resize(row, column);
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    void device_obj<DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>>::
            swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        Dim::swap(obj);
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    inline auto DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>::
            toDevice() const {
        return device_obj<This>(*this);
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    inline auto DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>::
            toDeviceAsync() const {
        device_obj<This> result(getRow(), getColumn());
        toDeviceAsync(result);
        return device_obj<This>(std::move(result));
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    inline void DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>::
            toDevice(device_obj<This>& obj) const {
        Base::toDevice(obj);
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    inline void DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>::
            toDeviceAsync(device_obj<This>& obj) const {
        Base::toDeviceAsync(obj);
    }
    //////////////////////////////////////////////Row-Element//////////////////////////////////////////////
    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    __host__ __device__ device_obj<DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>>::
            device_obj(size_t row, size_t column) {
        resize(row, column);
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    __host__ __device__ device_obj<DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>>::
            device_obj(size_t row, size_t column, T value) {
        resize(row, column, std::move(value));
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    device_obj<DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>>::
            device_obj(const host_obj& storage): Base(storage), Dim(storage.getRow(), storage.getColumn()) {}

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    __device__ T&
    device_obj<DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>>::
            operator()(size_t r, size_t c) {
        assert(r < getRow() && c < getColumn());
        return Base::operator[](toIndex(r, c));
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    __device__ const T&
    device_obj<DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>>::
            operator()(size_t r, size_t c) const {
        assert(r < getRow() && c < getColumn());
        return Base::operator[](toIndex(r, c));
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    template<class... Args>
    __host__ __device__ void device_obj<DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>>::
            resize(size_t row, size_t column, Args&&... args) {
    #ifdef __CUDA_ARCH__
        assert(Row * Column != Dynamic && "[Error]: Do not allocate dynamic matrix in device code");
    #endif
        Base::resize(row * column, std::forward<Args>(args)...);
        Dim::resize(row, column);
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    void device_obj<DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>>::
            swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        Dim::swap(obj);
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    inline auto DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>::
            toDevice() const {
        return device_obj<This>(*this);
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    inline auto DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>::
            toDeviceAsync() const {
        device_obj<This> result(getRow(), getColumn());
        toDeviceAsync(result);
        return device_obj<This>(std::move(result));
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    inline void DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>::
            toDevice(device_obj<This>& obj) const {
        Base::toDevice(obj);
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    inline void DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>::
            toDeviceAsync(device_obj<This>& obj) const {
        Base::toDeviceAsync(obj);
    }
    //////////////////////////////////////////////Column-Vector//////////////////////////////////////////////
    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    device_obj<DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>>::
            device_obj(size_t row, size_t column) : Dim(row, column), array(column, row), size(row * column) {}

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    device_obj<DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>>::
            device_obj(size_t row, size_t column, T value) : Dim(row, column), array(column, row, std::move(value)), size(row * column) {}

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    device_obj<DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>>::
            device_obj(const host_obj& storage) : Dim(storage.getRow(), storage.getColumn()), array(storage.array), size(storage.getSize()) {}

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    __device__ T&
    device_obj<DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>>::
            operator()(size_t r, size_t c) {
        assert(r < getRow() && c < getColumn());
        return array[c][r];
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    __device__ const T&
    device_obj<DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>>::
            operator()(size_t r, size_t c) const {
        assert(r < getRow() && c < getColumn());
        return array[c][r];
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    void device_obj<DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>>::
            resize(size_t row, size_t column) {
        array.resize(column);
        auto host_array = array.toHost();
        for (auto& vector : host_array)
            vector.resize(row);
        host_array.toDevice(array);
        Dim::resize(row, column);
        size = row * column;
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    inline void device_obj<DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>>::
            swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Dim::swap(obj);
        array.swap(obj.array);
        std::swap(size, obj.size);
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    __host__ __device__ T*
    device_obj<DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>>::
            data_ptr(size_t row, size_t column) {
        return const_cast<T*>(const_cast<const This&>(*this).data_ptr(row, column));
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    __host__ __device__ const T*
    device_obj<DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>>::
            data_ptr(size_t row, size_t column) const {
    #ifdef __CUDA_ARCH__
        return array[column].data() + row;
    #else
        const auto host_array = array.toPlainHost();
        const auto* p = host_array[column].getDerived().data() + row;
        return p;
    #endif
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    inline auto DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>::toDevice() const {
        return device_obj<This>(*this);
    }
    //////////////////////////////////////////////Row-Vector//////////////////////////////////////////////
    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    device_obj<DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>>::
            device_obj(size_t row, size_t column) : Dim(row, column), array(row, column), size(row * column) {}

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    device_obj<DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>>::
            device_obj(size_t row, size_t column, T value) : Dim(row, column), array(row, column, std::move(value)), size(row * column) {}

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    device_obj<DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>>::
            device_obj(const host_obj& storage) : Dim(storage.getRow(), storage.getColumn()), array(storage.array), size(storage.getSize()) {}

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    __device__ T&
    device_obj<DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>>::
            operator()(size_t r, size_t c) {
        assert(r < getRow() && c < getColumn());
        return array[r][c];
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    __device__ const T&
    device_obj<DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>>::
            operator()(size_t r, size_t c) const {
        assert(r < getRow() && c < getColumn());
        return array[r][c];
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    void device_obj<DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>>::
            resize(size_t row, size_t column) {
        array.resize(row);
        auto host_array = array.toHost();
        for (auto& vector : host_array)
            vector.resize(column);
        host_array.toDevice(array);
        Dim::resize(row, column);
        size = row * column;
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    inline void device_obj<DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>>::
            swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Dim::swap(obj);
        array.swap(obj.array);
        std::swap(size, obj.size);
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    __host__ __device__ T*
    device_obj<DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>>::
            data_ptr(size_t row, size_t column) {
        return const_cast<T*>(const_cast<const This&>(*this).data_ptr(row, column));
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    __host__ __device__ const T*
    device_obj<DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>>::
            data_ptr(size_t row, size_t column) const {
    #ifdef __CUDA_ARCH__
        return array[row].data() + column;
    #else
        const auto host_array = array.toPlainHost();
        const auto* p = host_array[row].getDerived().data() + column;
        return p;
    #endif
    }

    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    inline auto DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>::toDevice() const {
        return device_obj<This>(*this);
    }
}
