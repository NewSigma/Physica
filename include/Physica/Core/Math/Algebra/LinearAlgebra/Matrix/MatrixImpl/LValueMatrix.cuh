/*
 * Copyright 2022-2023 WeiBo He.
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

#include "RValueMatrix.cuh"
#include "LValueFlatten.cuh"

namespace Physica::Core {
    template<class Derived>
    class device_obj<LValueMatrix<Derived>> : public device_obj<RValueMatrix<Derived>> {
        using Base = device_obj<RValueMatrix<Derived>>;
    public:
        using typename Base::ScalarType;
        using Base::RowAtCompile;
        using Base::ColumnAtCompile;
    public:
        /* Operators */
        device_obj& operator=(const device_obj& m) = delete;
        device_obj& operator=(device_obj&& m) = delete;
        template<class OtherMatrix> __device__ device_obj<Derived>& operator=(const device_obj<RValueMatrix<OtherMatrix>>& m);
        __device__ Derived& operator=(const ScalarType& s);
        [[nodiscard]] __device__ ScalarType& operator()(size_t row, size_t column) { return Base::getDerived()(row, column); }
        [[nodiscard]] __device__ const ScalarType& operator()(size_t row, size_t column) const { return Base::getDerived()(row, column); }
        __device__ void operator+=(const ScalarType& s) { (*this) = (*this) + s; }
        __device__ void operator-=(const ScalarType& s) { (*this) = (*this) - s; }
        __device__ void operator*=(const ScalarType& s) { (*this) = (*this) * s; }
        __device__ void operator/=(const ScalarType& s) { (*this) = (*this) / s; }
        template<class OtherDerived> __device__ void operator+=(const device_obj<RValueMatrix<OtherDerived>>& m) { (*this) = (*this) + m; }
        template<class OtherDerived> __device__ void operator-=(const device_obj<RValueMatrix<OtherDerived>>& m) { (*this) = (*this) - m; }
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const { return (*this)(row, col); }
        [[nodiscard]] __host__ __device__ ScalarType* data_ptr(size_t row, size_t column) { return Base::getDerived().data_ptr(row, column); }
        [[nodiscard]] __host__ __device__ const ScalarType* data_ptr(size_t row, size_t column) const { return Base::getDerived().data_ptr(row, column); }
        [[nodiscard]] __device__ inline ScalarType& refFromMajorMinor(size_t major, size_t minor);
        [[nodiscard]] __device__ inline const ScalarType& refFromMajorMinor(size_t major, size_t minor) const;
        [[nodiscard]] __host__ __device__ device_obj<LValueFlatten<Derived>> flatten() { return {*this}; }
        [[nodiscard]] __host__ __device__ const device_obj<LValueFlatten<Derived>> flatten() const { return {*this}; }
    protected:
        device_obj() = default;
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
    };

    template<class Derived>
    template<class OtherMatrix>
    __device__ device_obj<Derived>& device_obj<LValueMatrix<Derived>>::operator=(const device_obj<RValueMatrix<OtherMatrix>>& m) {
        static_assert(RowAtCompile == Dynamic || OtherMatrix::RowAtCompile == Dynamic || RowAtCompile == OtherMatrix::RowAtCompile, "[Error]: Incompatible row number");
        static_assert(ColumnAtCompile == Dynamic || OtherMatrix::ColumnAtCompile == Dynamic || ColumnAtCompile == OtherMatrix::ColumnAtCompile, "[Error]: Incompatible column number");
        assert(getRow() == m.getRow() && "[Error]: Incompatible row number");
        assert(getColumn() == m.getColumn() && "[Error]: Incompatible column number");
        m.getDerived().assignTo(*this);
        return Base::getDerived();
    }

    template<class Derived>
    __device__ inline typename device_obj<LValueMatrix<Derived>>::ScalarType&
    device_obj<LValueMatrix<Derived>>::refFromMajorMinor(size_t major, size_t minor) {
        const size_t r = MatrixOption::rowFromMajorMinor<Derived>(major, minor);
        const size_t c = MatrixOption::columnFromMajorMinor<Derived>(major, minor);
        assert(r < Base::getDerived().getRow() && c < Base::getDerived().getColumn());
        return Base::getDerived()(r, c);
    }

    template<class Derived>
    __device__ inline const typename device_obj<LValueMatrix<Derived>>::ScalarType&
    device_obj<LValueMatrix<Derived>>::refFromMajorMinor(size_t major, size_t minor) const {
        const size_t r = MatrixOption::rowFromMajorMinor<Derived>(major, minor);
        const size_t c = MatrixOption::columnFromMajorMinor<Derived>(major, minor);
        assert(r < Base::getDerived().getRow() && c < Base::getDerived().getColumn());
        return Base::getDerived()(r, c);
    }
}
