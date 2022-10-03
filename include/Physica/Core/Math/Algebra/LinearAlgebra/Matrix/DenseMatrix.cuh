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

namespace Physica::Core {
    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    class device_obj<DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn>>
            : public device_obj<LValueMatrix<DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn>>>
            , public Internal::device_obj<Internal::DenseMatrixStorage<DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn>, option>>
            , public DenseMatrixDim<DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn>, Row, Column, MaxRow, MaxColumn> {
        using host_obj = DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn>;
        using This = device_obj<host_obj>;
        using Storage = Internal::device_obj<Internal::DenseMatrixStorage<host_obj, option>>;
        using Dim = DenseMatrixDim<host_obj, Row, Column, MaxRow, MaxColumn>;
    public:
        device_obj() = default;
        device_obj(const host_obj& mat);
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(device_obj obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] host_obj toHost() const { return host_obj(Storage::toHost()); }
        void swap(device_obj& obj) noexcept;
    };

    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    inline device_obj<DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn>>
    DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn>::toDevice() const {
        return device_obj<DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn>>(*this);
    }

    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    void device_obj<DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn>>::swap(
            device_obj<DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn>>& obj) noexcept {
        Storage::swap(obj);
        Dim::swap(obj);
    }
}
