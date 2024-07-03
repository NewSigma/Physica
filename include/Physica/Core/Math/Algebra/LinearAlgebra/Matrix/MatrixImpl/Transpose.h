/*
 * Copyright 2021-2024 WeiBo He.
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

#include "RValueMatrix.h"

namespace Physica::Core {
    template<class MatrixType>
    class Transpose : public RValueMatrix<Transpose<MatrixType>> {
        using Base = RValueMatrix<Transpose<MatrixType>>;

        const MatrixType& matrix;
    public:        
        using typename Base::ScalarType;
    public:
        Transpose(const RValueMatrix<MatrixType>& matrix_) : matrix(matrix_.getDerived()) {}
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return matrix.calc(col, row); }
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return matrix.getColumn(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const noexcept { return matrix.getRow(); }
    };

    template<class VectorType>
    class TransposeVector : public RValueMatrix<TransposeVector<VectorType>> {
        using This = TransposeVector<VectorType>;
    public:
        using Base = RValueMatrix<This>;
        using typename Base::ScalarType;
    private:
        const VectorType& vec;
    public:
        explicit TransposeVector(const RValueVector<VectorType>& vec_) : vec(vec_.getDerived()) {}
        /* Operations */
        template<class OtherMatrix>
        void assignTo(LValueMatrix<OtherMatrix>& target) const;
        /* Getters */
        [[nodiscard]] ScalarType calc([[maybe_unused]] size_t row, size_t col) const { assert(row == 0); return vec.calc(col); }
        [[nodiscard]] __host__ __device__ constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ size_t getColumn() const noexcept { return vec.getLength(); }
    };

    template<class VectorType>
    template<class OtherMatrix>
    void TransposeVector<VectorType>::assignTo(LValueMatrix<OtherMatrix>& target) const {
        using TargetType = LValueMatrix<OtherMatrix>;
        for (size_t i = 0; i < vec.getLength(); ++i)
            target.refFromMajorMinor(0, i) = calc(TargetType::rowFromMajorMinor(0, i), TargetType::columnFromMajorMinor(0, i));
    }
}

namespace Physica {
    using namespace Core;

    template<class MatrixType>
    class Traits<Transpose<MatrixType>> {
    private:
        constexpr static int OtherMajor = MatrixOption::isColumnMatrix<MatrixType>() ? MatrixOption::Row : MatrixOption::Column;
        constexpr static int Major = MatrixOption::isAnyMajor<MatrixType>() ? MatrixOption::AnyMajor : OtherMajor;
        constexpr static int Storage = MatrixOption::getStorage<MatrixType>();
    public:
        using ScalarType = typename MatrixType::ScalarType;
        constexpr static int Option = Major | Storage;
        constexpr static size_t RowAtCompile = MatrixType::ColumnAtCompile;
        constexpr static size_t ColumnAtCompile = MatrixType::RowAtCompile;
        constexpr static size_t MaxRowAtCompile = MatrixType::MaxColumnAtCompile;
        constexpr static size_t MaxColumnAtCompile = MatrixType::MaxRowAtCompile;
        constexpr static size_t SizeAtCompile = MatrixType::SizeAtCompile;
        constexpr static size_t MaxSizeAtCompile = MatrixType::MaxSizeAtCompile;
    };

    template<class VectorType>
    class Traits<TransposeVector<VectorType>> {
    public:
        using ScalarType = typename VectorType::ScalarType;
        constexpr static int Option = MatrixOption::Row | MatrixOption::Vector;
        constexpr static size_t RowAtCompile = 1;
        constexpr static size_t ColumnAtCompile = VectorType::SizeAtCompile;
        constexpr static size_t MaxRowAtCompile = 1;
        constexpr static size_t MaxColumnAtCompile = VectorType::MaxSizeAtCompile;
        constexpr static size_t SizeAtCompile = VectorType::SizeAtCompile;
        constexpr static size_t MaxSizeAtCompile = VectorType::MaxSizeAtCompile;
    };
}
