/*
 * Copyright 2022 WeiBo He.
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

#include "LValueMatrix.cuh"
#include "DenseMatrixImpl/DenseMatrixStorage.cuh"

namespace Physica::Core {
    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    class device_obj<DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>>
            : public device_obj<LValueMatrix<DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>>>
            , public Internal::device_obj<Internal::DenseMatrixStorage<DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>, option>>
            , public DenseMatrixDim<device_obj<DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>>, Row, Column, MaxRow, MaxColumn> {
        using host_obj = DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>;
        using host_storage = typename host_obj::Storage;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueMatrix<host_obj>>;
        using Storage = Internal::device_obj<Internal::DenseMatrixStorage<host_obj, option>>;
        using Dim = DenseMatrixDim<This, Row, Column, MaxRow, MaxColumn>;
    public:
        device_obj() = default;
        device_obj(const host_obj& mat);
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(device_obj obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] host_obj toHost() const { return host_obj(Storage::toHost(), Dim::getRow(), Dim::getColumn()); }
        void swap(device_obj& obj) noexcept;
        /* Getters */
        using Dim::getRow;
        using Dim::getColumn;
    };

    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    device_obj<DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>>::device_obj(const host_obj& mat)
            : Storage(static_cast<const host_storage&>(mat).toDevice()), Dim(mat.getRow(), mat.getColumn()) {}

    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    void device_obj<DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>>::swap(
            device_obj<DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>>& obj) noexcept {
        Storage::swap(obj);
        Dim::swap(obj);
    }

    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    inline device_obj<DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>>
    DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>::toDevice() const {
        return device_obj<DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn, Allocator>>(*this);
    }
}
