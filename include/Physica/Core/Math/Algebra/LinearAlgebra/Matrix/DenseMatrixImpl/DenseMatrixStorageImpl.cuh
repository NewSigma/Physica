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
    #define tparams class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator
    #define ColumnElementStorage DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>

    template<tparams>
    __host__ __device__ device_obj<ColumnElementStorage>::device_obj(size_t row, size_t column) {
        resize(row, column);
    }

    template<tparams>
    __host__ __device__ device_obj<ColumnElementStorage>::device_obj(size_t row, size_t column, ValueType value) {
        resize(row, column, std::move(value));
    }

    template<tparams>
    device_obj<ColumnElementStorage>::device_obj(const host_obj& storage) : Base(storage), Dim(storage.getRow(), storage.getColumn()) {}

    template<tparams>
    __device__ typename device_obj<ColumnElementStorage>::ValueType& device_obj<ColumnElementStorage>::operator()(size_t r, size_t c) {
        assert(r < getRow() && c < getColumn());
        return Base::operator[](toIndex(r, c));
    }

    template<tparams>
    __device__ const typename device_obj<ColumnElementStorage>::ValueType& device_obj<ColumnElementStorage>::operator()(size_t r, size_t c) const {
        assert(r < getRow() && c < getColumn());
        return Base::operator[](toIndex(r, c));
    }

    template<tparams>
    template<class... Args>
    __host__ __device__ void device_obj<ColumnElementStorage>::resize(size_t row, size_t column, Args&&... args) {
        if constexpr (IsDevice())
            assert(Row * Column != Dynamic && "[Error]: Do not allocate dynamic matrix in device code");
        Base::resize(row * column, std::forward<Args>(args)...);
        Dim::resize(row, column);
    }

    template<tparams>
    void device_obj<ColumnElementStorage>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        Dim::swap(obj);
    }

    template<tparams>
    inline auto ColumnElementStorage::toDevice() const {
        return device_obj<This>(*this);
    }

    template<tparams>
    inline auto ColumnElementStorage::toDeviceAsync() const {
        device_obj<This> result(getRow(), getColumn());
        toDeviceAsync(result);
        return device_obj<This>(std::move(result));
    }

    template<tparams>
    inline void ColumnElementStorage::toDevice(device_obj<This>& obj) const {
        Base::toDevice(obj);
    }

    template<tparams>
    inline void ColumnElementStorage::toDeviceAsync(device_obj<This>& obj) const {
        Base::toDeviceAsync(obj);
    }

    #undef ColumnElementStorage
    #define RowElementStorage DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>

    template<tparams>
    __host__ __device__ device_obj<RowElementStorage>::device_obj(size_t row, size_t column) {
        resize(row, column);
    }

    template<tparams>
    __host__ __device__ device_obj<RowElementStorage>::device_obj(size_t row, size_t column, ValueType value) {
        resize(row, column, std::move(value));
    }

    template<tparams>
    device_obj<RowElementStorage>::device_obj(const host_obj& storage): Base(storage), Dim(storage.getRow(), storage.getColumn()) {}

    template<tparams>
    __device__ typename device_obj<RowElementStorage>::ValueType&
    device_obj<RowElementStorage>::operator()(size_t r, size_t c) {
        assert(r < getRow() && c < getColumn());
        return Base::operator[](toIndex(r, c));
    }

    template<tparams>
    __device__ const typename device_obj<RowElementStorage>::ValueType&
    device_obj<RowElementStorage>::operator()(size_t r, size_t c) const {
        assert(r < getRow() && c < getColumn());
        return Base::operator[](toIndex(r, c));
    }

    template<tparams>
    template<class... Args>
    __host__ __device__ void device_obj<RowElementStorage>::resize(size_t row, size_t column, Args&&... args) {
        if constexpr (IsDevice())
            assert(Row * Column != Dynamic && "[Error]: Do not allocate dynamic matrix in device code");
        Base::resize(row * column, std::forward<Args>(args)...);
        Dim::resize(row, column);
    }

    template<tparams>
    void device_obj<RowElementStorage>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        Dim::swap(obj);
    }

    template<tparams>
    inline auto RowElementStorage::toDevice() const {
        return device_obj<This>(*this);
    }

    template<tparams>
    inline auto RowElementStorage::toDeviceAsync() const {
        device_obj<This> result(getRow(), getColumn());
        toDeviceAsync(result);
        return device_obj<This>(std::move(result));
    }

    template<tparams>
    inline void RowElementStorage::toDevice(device_obj<This>& obj) const {
        Base::toDevice(obj);
    }

    template<tparams>
    inline void RowElementStorage::toDeviceAsync(device_obj<This>& obj) const {
        Base::toDeviceAsync(obj);
    }

    #undef RowElementStorage
    #define ColumnVectorStorage DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>

    template<tparams>
    device_obj<ColumnVectorStorage>::device_obj(size_t row, size_t column)
            : Dim(row, column), array(column, row), size(row * column) {}

    template<tparams>
    device_obj<ColumnVectorStorage>::device_obj(size_t row, size_t column, T value)
            : Dim(row, column), array(column, row, std::move(value)), size(row * column) {}

    template<tparams>
    device_obj<ColumnVectorStorage>::device_obj(const host_obj& storage)
            : Dim(storage.getRow(), storage.getColumn()), array(storage.array), size(storage.getSize()) {}

    template<tparams>
    __device__ T&
    device_obj<ColumnVectorStorage>::operator()(size_t r, size_t c) {
        assert(r < getRow() && c < getColumn());
        return array[c][r];
    }

    template<tparams>
    __device__ const T&
    device_obj<ColumnVectorStorage>::operator()(size_t r, size_t c) const {
        assert(r < getRow() && c < getColumn());
        return array[c][r];
    }

    template<tparams>
    void device_obj<ColumnVectorStorage>::resize(size_t row, size_t column) {
        array.resize(column);
        auto host_array = array.toHost();
        for (auto& vector : host_array)
            vector.resize(row);
        host_array.toDevice(array);
        Dim::resize(row, column);
        size = row * column;
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
    device_obj<ColumnVectorStorage>::data_ptr(size_t row, size_t column) {
        return const_cast<T*>(const_cast<const This&>(*this).data_ptr(row, column));
    }

    template<tparams>
    __host__ __device__ const T*
    device_obj<ColumnVectorStorage>::data_ptr(size_t row, size_t column) const {
        if constexpr (IsDevice())
            return array[column].data() + row;
        else {
            const auto host_array = array.toPlainHost();
            const auto* p = host_array[column].getDerived().data() + row;
            return p;
        }
    }

    template<tparams>
    inline auto ColumnVectorStorage::toDevice() const {
        return device_obj<This>(*this);
    }

    #undef ColumnVectorStorage
    #define RowVectorStorage DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>

    template<tparams>
    device_obj<RowVectorStorage>::device_obj(size_t row, size_t column)
            : Dim(row, column), array(row, column), size(row * column) {}

    template<tparams>
    device_obj<RowVectorStorage>::device_obj(size_t row, size_t column, T value)
            : Dim(row, column), array(row, column, std::move(value)), size(row * column) {}

    template<tparams>
    device_obj<RowVectorStorage>::device_obj(const host_obj& storage)
            : Dim(storage.getRow(), storage.getColumn()), array(storage.array), size(storage.getSize()) {}

    template<tparams>
    __device__ T&
    device_obj<RowVectorStorage>::operator()(size_t r, size_t c) {
        assert(r < getRow() && c < getColumn());
        return array[r][c];
    }

    template<tparams>
    __device__ const T&
    device_obj<RowVectorStorage>::operator()(size_t r, size_t c) const {
        assert(r < getRow() && c < getColumn());
        return array[r][c];
    }

    template<tparams>
    void device_obj<RowVectorStorage>::resize(size_t row, size_t column) {
        array.resize(row);
        auto host_array = array.toHost();
        for (auto& vector : host_array)
            vector.resize(column);
        host_array.toDevice(array);
        Dim::resize(row, column);
        size = row * column;
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
    device_obj<RowVectorStorage>::data_ptr(size_t row, size_t column) {
        return const_cast<T*>(const_cast<const This&>(*this).data_ptr(row, column));
    }

    template<tparams>
    __host__ __device__ const T*
    device_obj<RowVectorStorage>::data_ptr(size_t row, size_t column) const {
        if constexpr (IsDevice())
            return array[row].data() + column;
        else {
            const auto host_array = array.toPlainHost();
            const auto* p = host_array[row].getDerived().data() + column;
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
