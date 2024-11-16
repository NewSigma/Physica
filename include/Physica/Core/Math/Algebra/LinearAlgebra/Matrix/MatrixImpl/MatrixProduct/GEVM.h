/*
 * Copyright 2021-2024 Weibo He.
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

#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/RValueMatrix.h>

namespace Physica::Core {
    template<class VectorType, class MatrixType>
    class VectorMatrixProduct : public RValueMatrix<VectorMatrixProduct<VectorType, MatrixType>> {
        using This = VectorMatrixProduct<VectorType, MatrixType>;
    public:
        using Base = RValueMatrix<This>;
        using Base::isReverseDiff;
        using typename Base::ScalarType;
    private:
        const VectorType& vec;
        const MatrixType& mat;
    public:
        VectorMatrixProduct(const RValueVector<VectorType>& vec_, const RValueMatrix<MatrixType>& mat_)
                : vec(vec_.getDerived()), mat(mat_.getDerived()) {
            assert(mat.getRow() == 1);
        }
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
        [[nodiscard]] const VectorType& getLHS() const noexcept { return vec; }
        [[nodiscard]] const MatrixType& getRHS() const noexcept { return mat; }
    };

    template<class VectorType, class MatrixType>
    typename VectorMatrixProduct<VectorType, MatrixType>::ScalarType VectorMatrixProduct<VectorType, MatrixType>::calc(size_t row, size_t col) const {
        return vec.calc(row) * mat.calc(0, col);
    }
    /**
     * \note Here we force the row of \param mat is 1, because in Physica vectors are naturally col vectors.
     * To compute row vector * matrix, users should converted it to matrix^T * col vector.
     */
    template<class VectorType, class MatrixType>
    [[nodiscard]] inline typename std::enable_if<MatrixType::RowAtCompile == 1, VectorMatrixProduct<VectorType, MatrixType>>::type
    operator*(const RValueVector<VectorType>& vec, const RValueMatrix<MatrixType>& mat) noexcept {
        assert(mat.getRow() == 1);
        return VectorMatrixProduct(vec, mat);
    }
}

namespace Physica {
    using namespace Core;

    template<class VectorType, class MatrixType>
    class Traits<VectorMatrixProduct<VectorType, MatrixType>> {
        static_assert(MatrixType::RowAtCompile == 1 || MatrixType::RowAtCompile == Dynamic,
                      "Row and column do not match in matrix product");
    public:
        using ScalarType = typename Core::Internal::BinaryScalarOpReturnType<typename VectorType::ScalarType,
                                                                             typename MatrixType::ScalarType>::Type;
        constexpr static int Option = MatrixOption::AnyMajor | MatrixOption::AnyStorage;
        constexpr static size_t RowAtCompile = VectorType::SizeAtCompile;
        constexpr static size_t ColAtCompile = MatrixType::ColAtCompile;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;
    };
}
