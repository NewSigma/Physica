/*
 * Copyright 2023 WeiBo He.
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
    namespace Internal {
        template<class T> class Traits;

        template<class MatrixType>
        class Traits<device_obj<ContinuousFlatten<MatrixType>>> : public Traits<ContinuousFlatten<MatrixType>> {};
    }

    template<class MatrixType>
    class device_obj<ContinuousFlatten<MatrixType>> : public device_obj<ContinuousVector<ContinuousFlatten<MatrixType>>> {
        using host_obj = ContinuousFlatten<MatrixType>;
        using This = device_obj<host_obj>;

        const device_obj<MatrixType>& mat;
    public:
        using Base = device_obj<ContinuousVector<host_obj>>;
        using typename Base::ScalarType;
    public:
        device_obj(const device_obj<ContinuousMatrix<MatrixType>>& mat_) : mat(mat_.getDerived()) {}
        /* Operators */
        [[nodiscard]] __device__ ScalarType& operator[](size_t index) { return *data_ptr(index); }
        [[nodiscard]] __device__ const ScalarType& operator[](size_t index) const { return *data_ptr(index); }
        /* Operations */
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return mat.getRow() * mat.getColumn(); }
        [[nodiscard]] __host__ __device__ inline ScalarType* data_ptr(size_t index);
        [[nodiscard]] __host__ __device__ inline const ScalarType* data_ptr(size_t index) const;
    };

    template<class MatrixType>
    __host__ __device__ inline typename device_obj<ContinuousFlatten<MatrixType>>::ScalarType*
    device_obj<ContinuousFlatten<MatrixType>>::data_ptr(size_t index) {
        return const_cast<ScalarType*>(const_cast<const This&>(*this).data_ptr(index));
    }

    template<class MatrixType>
    __host__ __device__ inline const typename device_obj<ContinuousFlatten<MatrixType>>::ScalarType*
    device_obj<ContinuousFlatten<MatrixType>>::data_ptr(size_t index) const {
        const size_t major = index / mat.getMaxMinor();
        const size_t minor = index % mat.getMaxMinor();
        const size_t row = MatrixOption::rowFromMajorMinor<MatrixType>(major, minor);
        const size_t col = MatrixOption::columnFromMajorMinor<MatrixType>(major, minor);
        return mat.data_ptr(row, col);
    }
}
