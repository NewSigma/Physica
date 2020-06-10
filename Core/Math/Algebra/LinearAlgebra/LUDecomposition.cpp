/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include "Core/Header/LUDecomposition.h"
#include "Core/Header/SquareMatrix.h"

namespace Physica::Core {
    LUDecomposition::LUDecomposition(SquareMatrix& square)
            : matrix(nullptr)
            , biasOrder(new size_t[square.getLength()]) {
        const auto length = square.getLength();
        for(size_t i = 0; i < length; ++i)
            biasOrder[i] = i;
        for(size_t i = 0; i < length; ++i) {
            std::swap(biasOrder[i], biasOrder[square.partialPivoting(i)]);
            decompositionColumn(square, square, i);
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
            decompositionColumn(*matrix, *matrix, i);
        }
    }

    LUDecomposition::~LUDecomposition() {
        delete matrix;
        delete[] biasOrder;
    }
    /*!
     * Apply LU Decomposition on a column of Matrix \from, save the result to Matrix \to.
     */
    void LUDecomposition::decompositionColumn(const Matrix& from, Matrix& to, size_t column) {
        const auto startAlphaIndex = column + 1;
        for (size_t j = 0; j < startAlphaIndex; ++j) {
            Numerical temp(from(j, column));
            for (size_t k = 0; k < j - 1; ++k)
                temp -= from(j, k) * from(k, column);
            to(j, column) = temp;
        }

        for (size_t j = startAlphaIndex; j < from.getLength(); ++j) {
            Numerical temp(from(j, column));
            for (size_t k = 0; k < j - 1; ++k)
                temp -= from(j, k) * from(k, column);
            to(j, column) = temp / from(column, column);
        }
    }
}