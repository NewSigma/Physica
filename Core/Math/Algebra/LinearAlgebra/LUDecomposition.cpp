/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include "Core/Header/LUDecomposition.h"
#include "Core/Header/SquareMatrix.h"

namespace Physica::Core {
    LUDecomposition::LUDecomposition(SquareMatrix& square)
            : matrix(nullptr)
            , biasOrder(new size_t[square.row()]) {
        const auto r = square.row();
        for(size_t i = 0; i < r; ++i)
            biasOrder[i] = i;
        for(size_t i = 0; i < r; ++i) {
            std::swap(biasOrder[i], biasOrder[square.partialPivoting(i)]);
            decompositionColumn(square, i);
        }
    }

    LUDecomposition::LUDecomposition(const SquareMatrix& square)
            : matrix(square.getType() == Matrix::Column
                     ? static_cast<Matrix*>(new ColumnSquareMatrix(static_cast<const ColumnSquareMatrix&>(square)))
                     : static_cast<Matrix*>(new RowSquareMatrix(static_cast<const RowSquareMatrix&>(square))))
            , biasOrder(new size_t[square.getLength()]) {
        const auto length = square.getLength();
        for (size_t i = 0; i < length; ++i)
            biasOrder[i] = i;
        for (size_t i = 0; i < length; ++i) {
            std::swap(biasOrder[i], biasOrder[matrix->partialPivoting(i)]);
            decompositionColumn(*matrix, i);
        }
    }

    LUDecomposition::~LUDecomposition() {
        delete matrix;
        delete[] biasOrder;
    }
    /*!
     * Apply LU Decomposition on a column of Matrix \from, save the result to Matrix \to.
     *
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009.32
     */
    void LUDecomposition::decompositionColumn(Matrix& m, size_t column) {
        const auto startAlphaIndex = column + 1;
        for (size_t j = 1; j < startAlphaIndex; ++j) {
            Scalar temp(m(j, column));
            for (size_t k = 0; k < j; ++k)
                temp -= m(j, k) * m(k, column);
            m(j, column) = std::move(temp);
        }

        const auto r = m.row();
        for (size_t j = startAlphaIndex; j < r; ++j) {
            Scalar temp(m(j, column));
            for (size_t k = 0; k < column; ++k)
                temp -= m(j, k) * m(k, column);
            m(j, column) = temp / m(column, column);
        }
    }
}