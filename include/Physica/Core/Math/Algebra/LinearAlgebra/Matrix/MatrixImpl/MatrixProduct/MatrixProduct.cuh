/*
 * Copyright 2023-2024 Weibo He.
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

#include <Physica/Utils/CUDA/PlainStruct.h>
#include "MatrixProduct.h"

namespace Physica::Core {
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
        Physica::PlainStruct<const device_obj<VectorType>> vec;
        Physica::PlainStruct<const device_obj<MatrixType>> mat;
    public:
        __host__ __device__ device_obj(const device_obj<RValueVector<VectorType>>& vec_, const device_obj<RValueMatrix<MatrixType>>& mat_)
                : vec(asStruct(vec_.getDerived())), mat(asStruct(mat_.getDerived())) {
            assert(mat.getDerived().getRow() == 1);
        }
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t column) const;
        [[nodiscard]] __host__ __device__ size_t getRow() const { return vec.getDerived().getLength(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return mat.getDerived().getColumn(); }
    };

    template<class VectorType, class MatrixType>
    __device__ typename device_obj<VectorMatrixProduct<VectorType, MatrixType>>::ScalarType
    device_obj<VectorMatrixProduct<VectorType, MatrixType>>::calc(size_t row, size_t column) const {
        return vec.getDerived().calc(row) * mat.getDerived().calc(0, column);
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
        Physica::PlainStruct<const device_obj<MatrixType>> mat;
        Physica::PlainStruct<const device_obj<VectorType>> vec;
    public:
        __host__ __device__ device_obj(const device_obj<RValueMatrix<MatrixType>>& mat_, const device_obj<RValueVector<VectorType>>& vec_)
                : mat(asStruct(mat_.getDerived())), vec(asStruct(vec_.getDerived())) {
            assert(mat.getDerived().getColumn() == vec.getDerived().getLength());
        }
        /* Getters */
        [[nodiscard]] __device__ inline ScalarType calc(size_t index) const;
        [[nodiscard]] __host__ __device__ size_t getLength() const { return mat.getDerived().getRow(); }
    };

    template<class MatrixType, class VectorType>
    __device__ inline typename device_obj<MatrixVectorProduct<MatrixType, VectorType>>::ScalarType
    device_obj<MatrixVectorProduct<MatrixType, VectorType>>::calc(size_t index) const {
        return mat.getDerived().row(index) * vec.getDerived();
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

namespace Physica {
    template<class VectorType, class MatrixType>
    class Traits<Core::device_obj<VectorMatrixProduct<VectorType, MatrixType>>>
            : public Traits<VectorMatrixProduct<VectorType, MatrixType>> {};

    template<class MatrixType, class VectorType>
    class Traits<Core::device_obj<MatrixVectorProduct<MatrixType, VectorType>>>
            : public Traits<MatrixVectorProduct<MatrixType, VectorType>> {};
}

#include "GEMM.cuh"
