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

#include <Physica/PlainStruct.h>

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
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const;
        [[nodiscard]] __host__ __device__ size_t getRow() const { return vec.getDerived().getLength(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const { return mat.getDerived().getCol(); }
    };

    template<class VectorType, class MatrixType>
    __device__ typename device_obj<VectorMatrixProduct<VectorType, MatrixType>>::ScalarType
    device_obj<VectorMatrixProduct<VectorType, MatrixType>>::calc(size_t row, size_t col) const {
        return vec.getDerived().calc(row) * mat.getDerived().calc(0, col);
    }

    template<class VectorType, class MatrixType>
    [[nodiscard]] __host__ __device__ inline typename std::enable_if<MatrixType::RowAtCompile == 1, device_obj<VectorMatrixProduct<VectorType, MatrixType>>>::type
    operator*(const device_obj<RValueVector<VectorType>>& vec, const device_obj<RValueMatrix<MatrixType>>& mat) noexcept {
        assert(mat.getRow() == 1);
        return {vec, mat};
    }
}

namespace Physica {
    template<class VectorType, class MatrixType>
    class Traits<Core::device_obj<VectorMatrixProduct<VectorType, MatrixType>>>
            : public Traits<VectorMatrixProduct<VectorType, MatrixType>> {};
}
