/*
 * Copyright 2020-2021 WeiBo He.
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
#include "DenseMatrixImpl/DenseMatrixBase.h"
#include "DenseMatrixImpl/DenseMatrixExpression.h"
#include "DenseMatrixImpl/MatrixProduct.h"
#include "MatrixDecomposition/Cholesky.h"
#include "MatrixDecomposition/LUDecomposition.h"
#include "InverseMatrix.h"

namespace Physica::Core {
    template<class T = MultiScalar, int type = DenseMatrixType::Column | DenseMatrixType::Vector
            , size_t Row = Dynamic, size_t Column = Dynamic, size_t MaxRow = Row, size_t MaxColumn = Column>
    class DenseMatrix;

    namespace Internal {
        template<class T>
        class Traits;

        template<class T, int type, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
        class Traits<DenseMatrix<T, type, Row, Column, MaxRow, MaxColumn>> {
        public:
            using ScalarType = T;
            constexpr static int MatrixType = type;
            constexpr static size_t RowAtCompile = Row;
            constexpr static size_t ColumnAtCompile = Column;
            constexpr static size_t MaxRowAtCompile = MaxRow;
            constexpr static size_t MaxColumnAtCompile = MaxColumn;
            constexpr static size_t SizeAtCompile = Row * Column;
            constexpr static size_t MaxSizeAtCompile = MaxRow * MaxColumn;
        };
    }

    template<class T, int type, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    std::ostream& operator<<(std::ostream& os, const DenseMatrix<T, type, Row, Column, MaxRow, MaxColumn>& mat);
    /**
     * DenseMatrix class
     * A matrix can be either fixed matrix, which have its max size defined,
     * or dynamic matrix, whose size is dynamically changed.
     * 
     * \tparam type
     * type is combinations of \enum DenseMatrixType
     */
    template<class T, int type, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    class DenseMatrix : public DenseMatrixBase<DenseMatrix<T, type, Row, Column, MaxRow, MaxColumn>> {
        using Base = DenseMatrixBase<DenseMatrix<T, type, Row, Column, MaxRow, MaxColumn>>;

        enum {
            isColumnMatrix = DenseMatrixType::isColumnMatrix<DenseMatrix>()
        };
    public:
        using Base::Base;
        template<Utils::ExpressionType expType, class T1, class T2>
        DenseMatrix(DenseMatrixExpression<expType, T1, T2> exp);
        template<class OtherMatrix>
        DenseMatrix(const DenseMatrixBase<OtherMatrix>& mat);
        template<class T1, class T2>
        DenseMatrix(MatrixProduct<T1, T2> pro);
        template<class MatrixIn>
        DenseMatrix(LUDecomposition<MatrixIn> lu);
        template<class MatrixIn>
        DenseMatrix(InverseMatrix<MatrixIn> inverse);
        /* Operators */
        using Base::operator=;
        template<Utils::ExpressionType expType, class T1, class T2>
        DenseMatrix& operator=(DenseMatrixExpression<expType, T1, T2> exp) { *this = std::move(DenseMatrix(exp)); return *this; }
        friend std::ostream& operator<<<>(std::ostream& os, const DenseMatrix& mat);
        /* Static members */
        [[nodiscard]] static DenseMatrix zeroMatrix(size_t row, size_t column) { return DenseMatrix(row, column, T(0)); }
        [[nodiscard]] static DenseMatrix unitMatrix(size_t order);
    };

    template<class T, int type, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    std::ostream& operator<<(std::ostream& os, const DenseMatrix<T, type, Row, Column, MaxRow, MaxColumn>& mat);

    template<class T, int type, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    inline void swap(DenseMatrix<T, type, Row, Column, MaxRow, MaxColumn>& m1
            , DenseMatrix<T, type, Row, Column, MaxRow, MaxColumn>& m2) noexcept { m1.swap(m2); }
}

#include "DenseMatrixImpl/DenseMatrixImpl.h"
