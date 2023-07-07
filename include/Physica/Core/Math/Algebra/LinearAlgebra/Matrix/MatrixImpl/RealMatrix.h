/*
 * Copyright 2022-2023 WeiBo He.
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
    template<class MatrixType> class RealMatrix;

    namespace Internal {
        template<class T> class Traits;

        template<class MatrixType>
        class Traits<RealMatrix<MatrixType>> {
        public:
            using ScalarType = typename MatrixType::ScalarType::RealType;
            constexpr static int Option = MatrixType::Option;
            constexpr static size_t RowAtCompile = MatrixType::RowAtCompile;
            constexpr static size_t ColumnAtCompile = MatrixType::ColumnAtCompile;
            constexpr static size_t MaxRowAtCompile = MatrixType::MaxRowAtCompile;
            constexpr static size_t MaxColumnAtCompile = MatrixType::MaxColumnAtCompile;
            constexpr static size_t SizeAtCompile = MatrixType::SizeAtCompile;
            constexpr static size_t MaxSizeAtCompile = MatrixType::MaxSizeAtCompile;
        };
    }

    template<class MatrixType>
    class RealMatrix : public RValueMatrix<RealMatrix<MatrixType>> {
    public:
        using ScalarType = typename MatrixType::ScalarType::RealType;
        using Base = RValueMatrix<RealMatrix<MatrixType>>;
    private:
        const MatrixType& mat;
    public:
        RealMatrix(const RValueMatrix<MatrixType>& mat_) : mat(mat_.getDerived()) {}

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return mat.calc(row, col).getReal(); }
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getColumn() const { return mat.getColumn(); }
    };

    template<class MatrixType>
    [[nodiscard]] RealMatrix<MatrixType> toRealMatrix(const RValueMatrix<MatrixType>& m) {
        return RealMatrix(m);
    }
}
