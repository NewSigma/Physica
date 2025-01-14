/*
 * Copyright 2021-2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/RValueMatrix.h"

namespace Physica::Core {
    template<Vector T, Matrix U>
    class VectorMatrixProduct : public RValueMatrix<VectorMatrixProduct<T, U>> {
        using This = VectorMatrixProduct<T, U>;
    public:
        using Base = RValueMatrix<This>;
        using Base::isReverseDiff;
        using typename Base::ScalarType;
    private:
        const T& vec;
        const U& mat;
    public:
        VectorMatrixProduct(const T& vec_, const U& mat_);
        VectorMatrixProduct(const This&) = delete;
        VectorMatrixProduct(This&&) noexcept = delete;
        ~VectorMatrixProduct() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const;
        [[nodiscard]] __host__ __device__ size_t getRow() const { return vec.getLength(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const { return mat.getCol(); }
        [[nodiscard]] const T& getLHS() const noexcept { return vec; }
        [[nodiscard]] const U& getRHS() const noexcept { return mat; }
    };

    template<Vector T, Matrix U>
    VectorMatrixProduct<T, U>::VectorMatrixProduct(const T& vec_, const U& mat_) : vec(vec_), mat(mat_) {
        assert(vec.getLength() > 0 && mat.getCol() > 0 && "[Error]: Empty vector or matrix");
        assert(mat.getRow() == 1 && "[Error]: Dimensions do not match");
    }

    template<Vector T, Matrix U>
    auto VectorMatrixProduct<T, U>::calc(size_t row, size_t col) const -> ScalarType {
        return vec.calc(row) * mat.calc(0, col);
    }
    /**
     * \note Here we force the row of \param mat is 1, because in Physica vectors are naturally col vectors.
     * To compute row vector * matrix, users should converted it to matrix^T * col vector.
     */
    template<Vector T, Matrix U>
    [[nodiscard]] inline auto operator*(const T& vec, const U& mat) noexcept requires(U::RowAtCompile == 1) {
        return VectorMatrixProduct<T, U>(vec, mat);
    }
}

namespace Physica {
    template<Vector T, Matrix U>
    class Traits<VectorMatrixProduct<T, U>> {
        static_assert(U::RowAtCompile == 1 || U::RowAtCompile == Dynamic,
                      "Row and column do not match in matrix product");
    public:
        using ScalarType = Core::Internal::BinaryScalarOpRtnTy<typename T::ScalarType, typename U::ScalarType>::Type;
        constexpr static int Option = MatrixOption::AnyMajor | MatrixOption::AnyStorage;
        constexpr static size_t RowAtCompile = T::SizeAtCompile;
        constexpr static size_t ColAtCompile = U::ColAtCompile;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;
    };
}
