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
    #define tparams class T, size_t Row, size_t Col, class Allocator
    #define ColumnElementStorage DenseMatrixStorage<T, MatrixOption::Col | MatrixOption::Element, Row, Col, Allocator>

    template<tparams>
    __host__ __device__ device_obj<ColumnElementStorage>::device_obj(size_t row, size_t col) : Dim(row, col), arr(row * col) {}

    template<tparams>
    __host__ __device__ device_obj<ColumnElementStorage>::device_obj(size_t row, size_t col, ElemType value) : Dim(row, col), arr(row * col, std::move(value)) {}

    template<tparams>
    device_obj<ColumnElementStorage>::device_obj(const host_obj& storage)
            : Dim(storage.getRow(), storage.getCol()), arr(storage.arr) {}

    template<tparams>
    __device__ device_obj<ColumnElementStorage>::ElemType& device_obj<ColumnElementStorage>::operator()(size_t r, size_t c) {
        assert(r < getRow() && c < getCol());
        return arr[toIndex(r, c)];
    }

    template<tparams>
    __device__ const device_obj<ColumnElementStorage>::ElemType& device_obj<ColumnElementStorage>::operator()(size_t r, size_t c) const {
        assert(r < getRow() && c < getCol());
        return arr[toIndex(r, c)];
    }

    template<tparams>
    template<class... Args>
    __host__ __device__ void device_obj<ColumnElementStorage>::resize(size_t row, size_t col, Args&&... args) {
        if constexpr (IsDevice())
            assert(Row * Col != Dynamic && "[Error]: Do not allocate dynamic matrix in device code");
        Dim::resize(row, col);
        arr.resize(row * col, std::forward<Args>(args)...);
    }

    template<tparams>
    void device_obj<ColumnElementStorage>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Dim::swap(obj);
        arr.swap(obj.arr);
    }

    template<tparams>
    inline auto ColumnElementStorage::toDevice() const {
        return device_obj<This>(*this);
    }

    template<tparams>
    inline auto ColumnElementStorage::toDeviceAsync() const {
        device_obj<This> result(getRow(), getCol());
        toDeviceAsync(result);
        return device_obj<This>(std::move(result));
    }

    template<tparams>
    inline void ColumnElementStorage::toDevice(device_obj<This>& obj) const {
        arr.toDevice(obj.asArray());
    }

    template<tparams>
    inline void ColumnElementStorage::toDeviceAsync(device_obj<This>& obj) const {
        arr.toDeviceAsync(obj.asArray());
    }

    #undef ColumnElementStorage
    #define RowElementStorage DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Element, Row, Col, Allocator>

    template<tparams>
    __host__ __device__ device_obj<RowElementStorage>::device_obj(size_t row, size_t col) : Dim(row, col), arr(row * col) {}

    template<tparams>
    __host__ __device__ device_obj<RowElementStorage>::device_obj(size_t row, size_t col, ElemType value) : Dim(row, col), arr(row * col, std::move(value)) {}

    template<tparams>
    device_obj<RowElementStorage>::device_obj(const host_obj& storage)
            : Dim(storage.getRow(), storage.getCol()), arr(storage.arr) {}

    template<tparams>
    __device__ device_obj<RowElementStorage>::ElemType&
    device_obj<RowElementStorage>::operator()(size_t r, size_t c) {
        assert(r < getRow() && c < getCol());
        return arr[toIndex(r, c)];
    }

    template<tparams>
    __device__ const device_obj<RowElementStorage>::ElemType&
    device_obj<RowElementStorage>::operator()(size_t r, size_t c) const {
        assert(r < getRow() && c < getCol());
        return arr[toIndex(r, c)];
    }

    template<tparams>
    template<class... Args>
    __host__ __device__ void device_obj<RowElementStorage>::resize(size_t row, size_t col, Args&&... args) {
        if constexpr (IsDevice())
            assert(Row * Col != Dynamic && "[Error]: Do not allocate dynamic matrix in device code");
        Dim::resize(row, col);
        arr.resize(row * col, std::forward<Args>(args)...);
    }

    template<tparams>
    void device_obj<RowElementStorage>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Dim::swap(obj);
        arr.swap(obj.arr);
    }

    template<tparams>
    inline auto RowElementStorage::toDevice() const {
        return device_obj<This>(*this);
    }

    template<tparams>
    inline auto RowElementStorage::toDeviceAsync() const {
        device_obj<This> result(getRow(), getCol());
        toDeviceAsync(result);
        return device_obj<This>(std::move(result));
    }

    template<tparams>
    inline void RowElementStorage::toDevice(device_obj<This>& obj) const {
        arr.toDevice(obj.asArray());
    }

    template<tparams>
    inline void RowElementStorage::toDeviceAsync(device_obj<This>& obj) const {
        arr.toDeviceAsync(obj.asArray());
    }

    #undef RowElementStorage
    #define ColumnVectorStorage DenseMatrixStorage<T, MatrixOption::Col | MatrixOption::Vector, Row, Col, Allocator>

    template<tparams>
    device_obj<ColumnVectorStorage>::device_obj(size_t row, size_t col)
            : Dim(row, col), array(col, row), size(row * col) {}

    template<tparams>
    device_obj<ColumnVectorStorage>::device_obj(size_t row, size_t col, T value)
            : Dim(row, col), array(col, row, std::move(value)), size(row * col) {}

    template<tparams>
    device_obj<ColumnVectorStorage>::device_obj(const host_obj& storage)
            : Dim(storage.getRow(), storage.getCol()), array(storage.array), size(storage.getSize()) {}

    template<tparams>
    __device__ T&
    device_obj<ColumnVectorStorage>::operator()(size_t r, size_t c) {
        assert(r < getRow() && c < getCol());
        return array[c][r];
    }

    template<tparams>
    __device__ const T&
    device_obj<ColumnVectorStorage>::operator()(size_t r, size_t c) const {
        assert(r < getRow() && c < getCol());
        return array[c][r];
    }

    template<tparams>
    void device_obj<ColumnVectorStorage>::resize(size_t row, size_t col) {
        array.resize(col);
        auto host_array = array.toHost();
        for (auto& vector : host_array)
            vector.resize(row);
        host_array.toDevice(array);
        Dim::resize(row, col);
        size = row * col;
    }

    template<tparams>
    inline void device_obj<ColumnVectorStorage>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Dim::swap(obj);
        array.swap(obj.array);
        std::swap(size, obj.size);
    }

    template<tparams>
    __host__ __device__ T*
    device_obj<ColumnVectorStorage>::data_ptr(size_t row, size_t col) {
        return const_cast<T*>(const_cast<const This&>(*this).data_ptr(row, col));
    }

    template<tparams>
    __host__ __device__ const T*
    device_obj<ColumnVectorStorage>::data_ptr(size_t row, size_t col) const {
        if constexpr (IsDevice())
            return array[col].data() + row;
        else {
            const auto host_array = array.toPlainHost();
            const auto* p = host_array[col].getDerived().data() + row;
            return p;
        }
    }

    template<tparams>
    inline auto ColumnVectorStorage::toDevice() const {
        return device_obj<This>(*this);
    }

    #undef ColumnVectorStorage
    #define RowVectorStorage DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Vector, Row, Col, Allocator>

    template<tparams>
    device_obj<RowVectorStorage>::device_obj(size_t row, size_t col)
            : Dim(row, col), array(row, col), size(row * col) {}

    template<tparams>
    device_obj<RowVectorStorage>::device_obj(size_t row, size_t col, T value)
            : Dim(row, col), array(row, col, std::move(value)), size(row * col) {}

    template<tparams>
    device_obj<RowVectorStorage>::device_obj(const host_obj& storage)
            : Dim(storage.getRow(), storage.getCol()), array(storage.array), size(storage.getSize()) {}

    template<tparams>
    __device__ T&
    device_obj<RowVectorStorage>::operator()(size_t r, size_t c) {
        assert(r < getRow() && c < getCol());
        return array[r][c];
    }

    template<tparams>
    __device__ const T&
    device_obj<RowVectorStorage>::operator()(size_t r, size_t c) const {
        assert(r < getRow() && c < getCol());
        return array[r][c];
    }

    template<tparams>
    void device_obj<RowVectorStorage>::resize(size_t row, size_t col) {
        array.resize(row);
        auto host_array = array.toHost();
        for (auto& vector : host_array)
            vector.resize(col);
        host_array.toDevice(array);
        Dim::resize(row, col);
        size = row * col;
    }

    template<tparams>
    inline void device_obj<RowVectorStorage>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Dim::swap(obj);
        array.swap(obj.array);
        std::swap(size, obj.size);
    }

    template<tparams>
    __host__ __device__ T*
    device_obj<RowVectorStorage>::data_ptr(size_t row, size_t col) {
        return const_cast<T*>(const_cast<const This&>(*this).data_ptr(row, col));
    }

    template<tparams>
    __host__ __device__ const T*
    device_obj<RowVectorStorage>::data_ptr(size_t row, size_t col) const {
        if constexpr (IsDevice())
            return array[row].data() + col;
        else {
            const auto host_array = array.toPlainHost();
            const auto* p = host_array[row].getDerived().data() + col;
            return p;
        }
    }

    template<tparams>
    inline auto RowVectorStorage::toDevice() const {
        return device_obj<This>(*this);
    }

    #undef RowVectorStorage
    #undef tparams
}
