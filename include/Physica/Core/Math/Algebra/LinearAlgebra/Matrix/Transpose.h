/*
 * Copyright 2021-2022 WeiBo He.
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
    template<class MatrixType> class Transpose;

    template<class VectorType> class TransposeVector;

    namespace Internal {
        template<class T> class Traits;

        template<class MatrixType>
        class Traits<Transpose<MatrixType>> {
        private:
            constexpr static int major = DenseMatrixOption::isColumnMatrix<MatrixType>() ? DenseMatrixOption::Row : DenseMatrixOption::Column;
            constexpr static int storage = DenseMatrixOption::getStorage<MatrixType>();
        public:
            using ScalarType = typename MatrixType::ScalarType;
            constexpr static int MatrixOption = major | storage;
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
            constexpr static int MatrixOption = DenseMatrixOption::Row | DenseMatrixOption::Vector;
            constexpr static size_t RowAtCompile = 1;
            constexpr static size_t ColumnAtCompile = VectorType::SizeAtCompile;
            constexpr static size_t MaxRowAtCompile = 1;
            constexpr static size_t MaxColumnAtCompile = VectorType::MaxSizeAtCompile;
            constexpr static size_t SizeAtCompile = VectorType::SizeAtCompile;
            constexpr static size_t MaxSizeAtCompile = VectorType::MaxSizeAtCompile;
        };
    }

    template<class MatrixType>
    class Transpose : public RValueMatrix<Transpose<MatrixType>> {
        const MatrixType& matrix;
    public:
        using Base = RValueMatrix<Transpose<MatrixType>>;
        using typename Base::ScalarType;
    public:
        Transpose(const RValueMatrix<MatrixType>& matrix_) : matrix(matrix_.getDerived()) {}
        template<class OtherMatrix>
        void assignTo(LValueMatrix<OtherMatrix>& target) const;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return matrix.calc(col, row); }
        [[nodiscard]] size_t getRow() const noexcept { return matrix.getColumn(); }
        [[nodiscard]] size_t getColumn() const noexcept { return matrix.getRow(); }
    };

    template<class MatrixType>
    template<class OtherMatrix>
    void Transpose<MatrixType>::assignTo(LValueMatrix<OtherMatrix>& target) const {
        using TargetType = LValueMatrix<OtherMatrix>;
        const size_t max_i = target.getMaxMajor();
        const size_t mat_j = target.getMaxMinor();
        for (size_t i = 0; i < max_i; ++i) {
            for (size_t j = 0; j < mat_j; ++j) {
                target.getElementFromMajorMinor(i, j) = calc(TargetType::rowFromMajorMinor(i, j),
                                                             TargetType::columnFromMajorMinor(i, j));
            }
        }
    }

    template<class VectorType>
    class TransposeVector : public RValueMatrix<TransposeVector<VectorType>> {
        const VectorType& vec;
    public:
        using Base = RValueMatrix<TransposeVector<VectorType>>;
        using typename Base::ScalarType;
    public:
        TransposeVector(const RValueVector<VectorType>& vec_) : vec(vec_.getDerived()) {}
        template<class OtherMatrix>
        void assignTo(LValueMatrix<OtherMatrix>& target) const;
        /* Getters */
        [[nodiscard]] ScalarType calc([[maybe_unused]] size_t row, size_t col) const { assert(row == 0); return vec.calc(col); }
        [[nodiscard]] constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] size_t getColumn() const noexcept { return vec.getLength(); }
    };

    template<class VectorType>
    template<class OtherMatrix>
    void TransposeVector<VectorType>::assignTo(LValueMatrix<OtherMatrix>& target) const {
        using TargetType = LValueMatrix<OtherMatrix>;
        for (size_t i = 0; i < vec.getLength(); ++i) {
            target.getElementFromMajorMinor(0, i) = calc(TargetType::rowFromMajorMinor(0, i),
                                                         TargetType::columnFromMajorMinor(0, i));
        }
    }
}
