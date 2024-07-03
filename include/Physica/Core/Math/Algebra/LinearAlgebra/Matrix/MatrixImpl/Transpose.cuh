/*
 * Copyright 2023-2024 WeiBo He.
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

#include "Transpose.h"
#include "RValueMatrix.cuh"

namespace Physica::Core {
    template<class MatrixType>
    class device_obj<Transpose<MatrixType>> : public device_obj<RValueMatrix<Transpose<MatrixType>>> {
        using host_obj = Transpose<MatrixType>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;

        const device_obj<MatrixType>& mat;
    public:
        using typename Base::ScalarType;
    public:
        __host__ __device__ device_obj(const device_obj<RValueMatrix<MatrixType>>& mat_) : mat(mat_.getDerived()) {}
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const { return mat.calc(col, row); }
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return mat.getColumn(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const noexcept { return mat.getRow(); }
    };

    template<class VectorType>
    class device_obj<TransposeVector<VectorType>> : public device_obj<RValueMatrix<TransposeVector<VectorType>>> {
        using host_obj = TransposeVector<VectorType>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;

        const device_obj<VectorType>& vec;
    public:
        using typename Base::ScalarType;
    public:
        __host__ __device__ explicit device_obj(const device_obj<RValueVector<VectorType>>& vec_) : vec(vec_.getDerived()) {}
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc([[maybe_unused]] size_t row, size_t col) const { assert(row == 0); return vec.calc(col); }
        [[nodiscard]] __host__ __device__ constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ size_t getColumn() const noexcept { return vec.getLength(); }
    };
}

namespace Physica {
    using namespace Core;

    template<class MatrixType>
    class Traits<Core::device_obj<Transpose<MatrixType>>> : public Traits<Transpose<MatrixType>> {};

    template<class VectorType>
    class Traits<Core::device_obj<TransposeVector<VectorType>>> : public Traits<TransposeVector<VectorType>> {};
}
