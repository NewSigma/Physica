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
        [[nodiscard]] ScalarType calc(size_t row, size_t column) const;
        [[nodiscard]] __host__ __device__ size_t getRow() const { return vec.getLength(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return mat.getColumn(); }
        [[nodiscard]] const VectorType& getLHS() const noexcept { return vec; }
        [[nodiscard]] const MatrixType& getRHS() const noexcept { return mat; }
    };

    template<class MatrixType, class VectorType>
    class MatrixVectorProduct : public RValueVector<MatrixVectorProduct<MatrixType, VectorType>> {
        using This = MatrixVectorProduct<MatrixType, VectorType>;
    public:
        using Base = RValueVector<This>;
        using typename Base::ScalarType;
    private:
        const MatrixType& mat;
        const VectorType& vec;
    public:
        MatrixVectorProduct(const RValueMatrix<MatrixType>& mat_, const RValueVector<VectorType>& vec_)
                : mat(mat_.getDerived()), vec(vec_.getDerived()) {
            assert(mat.getColumn() == vec.getLength());
        }
        MatrixVectorProduct(const This&) = delete;
        MatrixVectorProduct(This&&) noexcept = delete;
        ~MatrixVectorProduct() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<class OtherDerived, class Executor = SequentialExecutor>
        inline void assignTo(LValueVector<OtherDerived>& target) const;
        /* Getters */
        [[nodiscard]] inline ScalarType calc(size_t index) const;
        [[nodiscard]] __host__ __device__ size_t getLength() const { return mat.getRow(); }
        [[nodiscard]] const MatrixType& getLHS() const noexcept { return mat; }
        [[nodiscard]] const VectorType& getRHS() const noexcept { return vec; }
    };

    template<class VectorType, class MatrixType>
    typename VectorMatrixProduct<VectorType, MatrixType>::ScalarType VectorMatrixProduct<VectorType, MatrixType>::calc(size_t row, size_t column) const {
        return vec.calc(row) * mat.calc(0, column);
    }

    template<class MatrixType, class VectorType>
    template<class OtherDerived, class Executor>
    inline void MatrixVectorProduct<MatrixType, VectorType>::assignTo(LValueVector<OtherDerived>& target) const {
        for (size_t i = 0; i < getLength(); ++i)
            target[i] = calc(i);
        
        constexpr bool isContinuous = std::is_base_of<ContinuousVector<OtherDerived>, OtherDerived>::value;
        if constexpr (isContinuous && Base::isReverseDiff)
            target.getDerived().makeContinuous();
    }

    template<class MatrixType, class VectorType>
    inline typename MatrixVectorProduct<MatrixType, VectorType>::ScalarType MatrixVectorProduct<MatrixType, VectorType>::calc(size_t index) const {
        return mat.row(index) * vec;
    }

    /**
     * \note Here we force the row of \param mat is 1, because in Physica vectors are naturally column vectors.
     * To compute row vector * matrix, users should converted it to matrix^T * column vector.
     */
    template<class VectorType, class MatrixType>
    [[nodiscard]] inline typename std::enable_if<MatrixType::RowAtCompile == 1, VectorMatrixProduct<VectorType, MatrixType>>::type
    operator*(const RValueVector<VectorType>& vec, const RValueMatrix<MatrixType>& mat) noexcept {
        assert(mat.getRow() == 1);
        return VectorMatrixProduct(vec, mat);
    }

    template<class MatrixType, class VectorType>
    [[nodiscard]] inline typename std::enable_if<MatrixType::RowAtCompile != 1, MatrixVectorProduct<MatrixType, VectorType>>::type
    operator*(const RValueMatrix<MatrixType>& mat, const RValueVector<VectorType>& vec) noexcept {
        assert(mat.getColumn() == vec.getLength());
        return MatrixVectorProduct(mat.getDerived(), vec.getDerived());
    }

    template<class MatrixType, class VectorType>
    [[nodiscard]] inline typename std::enable_if<MatrixType::RowAtCompile == 1 && MatrixType::ColumnAtCompile == 1,
                                                 typename Internal::BinaryScalarOpReturnType<typename MatrixType::ScalarType,
                                                                                             typename VectorType::ScalarType>::Type>::type
    operator*(const RValueMatrix<MatrixType>& mat, const RValueVector<VectorType>& vec) {
        assert(mat.getColumn() == vec.getLength());
        return mat.row(0) * vec;
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
        constexpr static size_t ColumnAtCompile = MatrixType::ColumnAtCompile;
        constexpr static size_t MaxRowAtCompile = VectorType::MaxSizeAtCompile;
        constexpr static size_t MaxColumnAtCompile = MatrixType::MaxColumnAtCompile;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColumnAtCompile;
        constexpr static size_t MaxSizeAtCompile = SizeAtCompile;
    };

    template<class MatrixType, class VectorType>
    class Traits<MatrixVectorProduct<MatrixType, VectorType>> {
        static_assert(MatrixType::ColumnAtCompile == VectorType::SizeAtCompile ||
                      MatrixType::ColumnAtCompile == Dynamic ||
                      VectorType::SizeAtCompile == Dynamic,
                      "Row and column do not match in matrix product");
    public:
        using ScalarType = typename Core::Internal::BinaryScalarOpReturnType<typename MatrixType::ScalarType,
                                                                             typename VectorType::ScalarType>::Type;
        constexpr static size_t SizeAtCompile = MatrixType::RowAtCompile;
        constexpr static size_t MaxSizeAtCompile = MatrixType::MaxRowAtCompile;

        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = false;
    };
}

#include "GEMM.h"
