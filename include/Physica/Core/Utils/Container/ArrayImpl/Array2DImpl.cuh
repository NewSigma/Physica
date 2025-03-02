/*
 * Copyright 2024-2025 Weibo He.
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

#include "../Array2D.cuh"

namespace Physica {
#define tparams class T, int Option, size_t Row, size_t Col, class Allocator
#define Array2D Array2D<T, Option, Row, Col, Allocator>

    template<tparams>
    __host__ __device__ device_obj<Array2D>::device_obj(size_t row, size_t col) {
        resize(row, col);
    }

    template<tparams>
    __host__ __device__ device_obj<Array2D>::device_obj(size_t row, size_t col, ElemType value) {
        resize(row, col, std::move(value));
    }

    template<tparams>
    device_obj<Array2D>::device_obj(const host_obj& storage)
            : arr(storage.arr), r(storage.r), size(storage.getSize()) {}

    template<tparams>
    __device__ auto device_obj<Array2D>::operator()(size_t r, size_t c) -> ElemType& {
        assert(r < getRow() && c < getCol());
        return *data_ptr(r, c);
    }

    template<tparams>
    __device__ auto device_obj<Array2D>::operator()(size_t r, size_t c) const -> const ElemType& {
        return const_cast<This&>(*this).operator()(r, c);
    }

    template<tparams>
    template<class... Args>
    __host__ __device__ void device_obj<Array2D>::resize(size_t row, size_t col, Args&&... args) {
        assert((Row == row || Row == Dynamic) && "[Error]: Cannot resize a fixed array");
        assert((Col == col || Col == Dynamic) && "[Error]: Cannot resize a fixed array");
        assert(row > 0 && col > 0);
        if constexpr (isVectorStorage) {
            if constexpr (IsDevice())
                noImpl("[Error]: resize() is not supported on device");
            else {
                const size_t major = isColMajor ? col : row;
                const size_t minor = isColMajor ? row : col;
                arr.resize(major);
                auto host_arr = arr.toHost();
                for (auto& elem : host_arr)
                    elem.resize(minor, std::forward<Args>(args)...);
                host_arr.toDevice(arr);
                size = row * col;
            }
        }
        else {
            if constexpr (IsDevice())
                assert(Row * Col != Dynamic && "[Error]: Do not allocate dynamic matrix in device code");
            else
                arr.resize(row * col, std::forward<Args>(args)...);
        }
        r = row;
    }

    template<tparams>
    auto device_obj<Array2D>::toHost() const -> host_obj {
        return host_obj(arr.toHost(), r);
    }

    template<tparams>
    void device_obj<Array2D>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        arr.swap(obj.arr);
        std::swap(r, obj.r);
        std::swap(size, obj.size);
    }

    template<tparams>
    __host__ __device__ size_t device_obj<Array2D>::getRow() const noexcept {
        if constexpr (Row == Dynamic) {
            if constexpr (Col == Dynamic)
                return r;
            else {
                const size_t size = getSize();
                assert(size % Col == 0);
                return size / Col;
            }
        }
        else
            return Row;
    }

    template<tparams>
    __host__ __device__ size_t device_obj<Array2D>::getCol() const noexcept {
        if constexpr (Col == Dynamic) {
            const size_t size = getSize();
            if constexpr (Row == Dynamic) {
                assert(r == 0 || size % getRow() == 0);
                return r == 0 ? 0 : size / getRow();
            }
            else {
                assert(size % Row == 0);
                return size / Row;
            }
        }
        else
            return Col;
    }

    template<tparams>
    __host__ __device__ size_t device_obj<Array2D>::getSize() const noexcept {
        if constexpr (isVectorStorage)
            return size;
        else
            return arr.getLength();
    }

    template<tparams>
    __host__ __device__ auto device_obj<Array2D>::data_ptr(size_t row, size_t col) -> ElemType* {
        assert(row < getRow());
        assert(col < getCol());
        if constexpr (isVectorStorage) {
            const size_t major = isColMajor ? col : row;
            const size_t minor = isColMajor ? row : col;
            if constexpr (IsDevice())
                return arr[major].data() + minor;
            else {
                auto host_arr = arr.toPlainHost();
                return host_arr[major].getDerived().data() + minor;
            }
        }
        else
            return arr.data() + toIndex1D(row, col);
    }

    template<tparams>
    __host__ __device__ auto device_obj<Array2D>::data_ptr(size_t row, size_t col) const -> const ElemType* {
        return const_cast<This&>(*this).data_ptr(row, col);
    }

    template<tparams>
    __host__ __device__ size_t device_obj<Array2D>::toIndex1D(size_t r, size_t c) const {
        if constexpr (isColMajor)
            return getRow() * c + r;
        else
            return getCol() * r + c;
    }

    template<tparams>
    inline auto Array2D::toDevice() const {
        return device_obj<This>(*this);
    }

    template<tparams>
    inline auto Array2D::toDeviceAsync() const {
        device_obj<This> result(getRow(), getCol());
        toDeviceAsync(result);
        return device_obj<This>(std::move(result));
    }

    template<tparams>
    inline void Array2D::toDevice(device_obj<This>& obj) const {
        arr.toDevice(obj.asArray());
    }

    template<tparams>
    inline void Array2D::toDeviceAsync(device_obj<This>& obj) const {
        arr.toDeviceAsync(obj.asArray());
    }

#undef Array2D
#undef tparams
}
