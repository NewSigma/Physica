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

#include "MatrixProduct.h"

namespace Physica::Core {
    namespace Internal {
        template<class VectorType, class MatrixType>
        class Traits<Core::device_obj<VectorMatrixProduct<VectorType, MatrixType>>>
                : public Traits<VectorMatrixProduct<VectorType, MatrixType>> {};

        template<class MatrixType, class VectorType>
        class Traits<Core::device_obj<MatrixVectorProduct<MatrixType, VectorType>>>
                : public Traits<MatrixVectorProduct<MatrixType, VectorType>> {};
    }

    template<class VectorType, class MatrixType>
    class device_obj<VectorMatrixProduct<VectorType, MatrixType>>
            : public device_obj<RValueMatrix<VectorMatrixProduct<VectorType, MatrixType>>> {
        static_assert(MatrixType::RowAtCompile == 1 || MatrixType::RowAtCompile == Dynamic,
                      "Row and column do not match in matrix product");
        using host_obj = VectorMatrixProduct<VectorType, MatrixType>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
    public:
        using typename Base::ScalarType;
    private:
        const device_obj<VectorType>& vec;
        const device_obj<MatrixType>& mat;
    public:
        __host__ __device__ device_obj(const device_obj<RValueVector<VectorType>>& vec_, const device_obj<RValueMatrix<MatrixType>>& mat_)
                : vec(vec_.getDerived()), mat(mat_.getDerived()) {
            assert(mat.getRow() == 1);
        }
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t column) const;
        [[nodiscard]] __host__ __device__ size_t getRow() const { return vec.getLength(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return mat.getColumn(); }
        [[nodiscard]] const VectorType& getLHS() const noexcept { return vec; }
        [[nodiscard]] const MatrixType& getRHS() const noexcept { return mat; }
    };

    template<class VectorType, class MatrixType>
    __device__ typename device_obj<VectorMatrixProduct<VectorType, MatrixType>>::ScalarType
    device_obj<VectorMatrixProduct<VectorType, MatrixType>>::calc(size_t row, size_t column) const {
        return vec.calc(row) * mat.calc(0, column);
    }

    template<class MatrixType, class VectorType>
    class device_obj<MatrixVectorProduct<MatrixType, VectorType>>
            : public device_obj<RValueVector<MatrixVectorProduct<MatrixType, VectorType>>> {
        static_assert(MatrixType::ColumnAtCompile == VectorType::SizeAtCompile ||
                      MatrixType::ColumnAtCompile == Dynamic ||
                      VectorType::SizeAtCompile == Dynamic,
                      "[Error]: Row and column do not match in matrix-vector product");
    public:
        using host_obj = MatrixVectorProduct<MatrixType, VectorType>;
        using Base = device_obj<RValueVector<host_obj>>;
        using typename Base::ScalarType;
    private:
        using This = device_obj<host_obj>;
        const device_obj<MatrixType>& mat;
        const device_obj<VectorType>& vec;
    public:
        __host__ __device__ device_obj(const device_obj<RValueMatrix<MatrixType>>& mat_, const device_obj<RValueVector<VectorType>>& vec_)
                : mat(mat_.getDerived()), vec(vec_.getDerived()) {
            assert(mat.getColumn() == vec.getLength());
        }
        /* Operations */
        template<class OtherDerived>
        __host__ __device__ inline void assignTo(device_obj<LValueVector<OtherDerived>>& target) const;
        /* Getters */
        [[nodiscard]] __device__ inline ScalarType calc(size_t index) const;
        [[nodiscard]] __host__ __device__ size_t getLength() const { return mat.getRow(); }
        [[nodiscard]] __host__ __device__ const MatrixType& getLHS() const noexcept { return mat; }
        [[nodiscard]] __host__ __device__ const VectorType& getRHS() const noexcept { return vec; }
    };

    template<class MatrixType, class VectorType>
    template<class OtherDerived>
    __host__ __device__ inline void device_obj<MatrixVectorProduct<MatrixType, VectorType>>::assignTo(device_obj<LValueVector<OtherDerived>>& target) const {
        for (size_t i = 0; i < getLength(); ++i)
            target[i] = calc(i);
    }

    template<class MatrixType, class VectorType>
    __device__ inline typename device_obj<MatrixVectorProduct<MatrixType, VectorType>>::ScalarType
    device_obj<MatrixVectorProduct<MatrixType, VectorType>>::calc(size_t index) const {
        return mat.row(index) * vec;
    }

    template<class VectorType, class MatrixType>
    __host__ __device__ inline typename std::enable_if<MatrixType::RowAtCompile == 1, device_obj<VectorMatrixProduct<VectorType, MatrixType>>>::type
    operator*(const device_obj<RValueVector<VectorType>>& vec, const device_obj<RValueMatrix<MatrixType>>& mat) {
        assert(mat.getRow() == 1);
        return {vec, mat};
    }

    template<class MatrixType, class VectorType>
    __host__ __device__ inline typename std::enable_if<MatrixType::RowAtCompile != 1, device_obj<MatrixVectorProduct<MatrixType, VectorType>>>::type
    operator*(const device_obj<RValueMatrix<MatrixType>>& mat, const device_obj<RValueVector<VectorType>>& vec) {
        assert(mat.getColumn() == vec.getLength());
        return {mat.getDerived(), vec.getDerived()};
    }
}
