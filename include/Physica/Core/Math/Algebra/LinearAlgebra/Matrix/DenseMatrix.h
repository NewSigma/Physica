/*
 * Copyright 2020-2022 WeiBo He.
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

#include <memory>
#include "LValueMatrix.h"
#include "DenseMatrixImpl/DenseMatrixExpression.h"
#include "DenseMatrixImpl/MatrixProduct.h"
#include "InverseMatrix.h"
#include "Transpose.h"
#include "MatrixDecomposition/LUDecomposition.h"

namespace Physica::Core {
    template<class T = MultiScalar, int option = DenseMatrixOption::Column | DenseMatrixOption::Vector
            , size_t Row = Dynamic, size_t Column = Dynamic, size_t MaxRow = Row, size_t MaxColumn = Column>
    class DenseMatrix;

    namespace Internal {
        template<class T>
        class Traits;

        template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
        class Traits<DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn>> {
        public:
            using ScalarType = T;
            constexpr static int MatrixOption = option;
            constexpr static size_t RowAtCompile = Row;
            constexpr static size_t ColumnAtCompile = Column;
            constexpr static size_t MaxRowAtCompile = MaxRow;
            constexpr static size_t MaxColumnAtCompile = MaxColumn;
            constexpr static size_t SizeAtCompile = Row * Column;
            constexpr static size_t MaxSizeAtCompile = MaxRow * MaxColumn;
        };
    }

    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    std::ostream& operator<<(std::ostream& os, const DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn>& mat);
    /**
     * DenseMatrix class
     * A matrix can be either fixed matrix, which have its max size defined,
     * or dynamic matrix, whose size is dynamically changed.
     * 
     * \tparam option
     * option is combinations of \enum DenseMatrixOption
     */
    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    class DenseMatrix : public LValueMatrix<DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn>>,
                        public DenseMatrixStorage<T, option, Row, Column, MaxRow, MaxColumn> {
        using Base = LValueMatrix<DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn>>;
        using Storage = DenseMatrixStorage<T, option, Row, Column, MaxRow, MaxColumn>;
        static_assert(MaxRow * MaxColumn * sizeof(T) <= 2048, "[Warning]: It is suggested declare large fixed size matrix as dynamic matrix");
    public:
        using ColMatrix = DenseMatrix<T, DenseMatrixOption::getStorage<DenseMatrix>() | DenseMatrixOption::Column, Row, Column, MaxRow, MaxColumn>;
        using RowMatrix = DenseMatrix<T, DenseMatrixOption::getStorage<DenseMatrix>() | DenseMatrixOption::Row, Row, Column, MaxRow, MaxColumn>;
        using RealMatrix = DenseMatrix<typename T::RealType, option, Row, Column, MaxRow, MaxColumn>;
    public:
        using Storage::Storage;
        template<class OtherMatrix>
        DenseMatrix(const RValueMatrix<OtherMatrix>& mat);
        template<class VectorType>
        DenseMatrix(const RValueVector<VectorType>& mat);
        template<class MatrixIn>
        DenseMatrix(LUDecomposition<MatrixIn> lu);
        DenseMatrix(const DenseMatrix& m);
        DenseMatrix(DenseMatrix&& m) noexcept;
        /* Operators */
        DenseMatrix& operator=(DenseMatrix m) noexcept;
        using Base::operator=;
        using Storage::operator();
        friend std::ostream& operator<<<>(std::ostream& os, const DenseMatrix& mat);
        /* Getters */
        using Storage::getRow;
        using Storage::getColumn;
        /* Helpers */
        void swap(DenseMatrix& m) noexcept;
        template<class VectorType>
        [[nodiscard]] static std::pair<DenseMatrix, DenseMatrix> meshgrid(const LValueVector<VectorType>& vecInCols, const LValueVector<VectorType>& vecInRows);
        /* Static members */
        [[nodiscard]] static DenseMatrix Zeros(size_t rank) { return DenseMatrix(rank, rank, T(0)); }
        [[nodiscard]] static DenseMatrix Zeros(size_t row, size_t column) { return DenseMatrix(row, column, T(0)); }
        [[nodiscard]] static DenseMatrix unitMatrix(size_t order);
        [[nodiscard]] static DenseMatrix randomMatrix(size_t order) { return randomMatrix(order, order); }
        [[nodiscard]] static DenseMatrix randomMatrix(size_t row, size_t column);
    };

    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    std::ostream& operator<<(std::ostream& os, const DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn>& mat);

    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    std::istream& operator>>(std::istream& is, DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn>& mat);

    template<class T, int option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    inline void swap(
            DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn>& m1,
            DenseMatrix<T, option, Row, Column, MaxRow, MaxColumn>& m2) noexcept {
        m1.swap(m2);
    }
}

#include "DenseMatrixImpl/DenseMatrixImpl.h"
