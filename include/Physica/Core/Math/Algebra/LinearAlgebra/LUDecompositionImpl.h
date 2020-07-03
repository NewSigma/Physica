/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_LUDECOMPOSITIONIMPL_H
#define PHYSICA_LUDECOMPOSITIONIMPL_H

#include "Matrix/MatrixOperation.h"

namespace Physica::Core {
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    LUDecomposition<T, type, maxRow, maxColumn>::LUDecomposition(const Matrix<T, type, maxRow, maxColumn>& m)
            : matrix(m), biasOrder(reinterpret_cast<size_t*>(malloc(m.getRow() * sizeof(size_t)))) {
        const auto rank = matrix.getRow();
        for(size_t i = 0; i < rank; ++i)
            biasOrder[i] = i;
        for(size_t i = 0; i < rank; ++i) {
            std::swap(biasOrder[i], biasOrder[MatrixOperation<T, type, maxRow, maxColumn>::partialPivoting(matrix, i)]);
            decompositionColumn(matrix, i);
        }
    }

    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    LUDecomposition<T, type, maxRow, maxColumn>::LUDecomposition(Matrix<T, type, maxRow, maxColumn>&& m) noexcept
            : matrix(std::move(m)), biasOrder(reinterpret_cast<size_t*>(malloc(m.getRow() * sizeof(size_t)))) {}

    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    LUDecomposition<T, type, maxRow, maxColumn>::LUDecomposition(const LUDecomposition& l)
            : matrix(l.matrix), biasOrder(reinterpret_cast<size_t*>(malloc(l.matrix.getRow() * sizeof(size_t)))) {
        const auto rank = l.matrix.getRow();
        for(size_t i = 0; i < rank; ++i) {
            biasOrder[i] = l.biasOrder[i];
        }
    }

    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    LUDecomposition<T, type, maxRow, maxColumn>::LUDecomposition(LUDecomposition&& l) noexcept
            : matrix(std::move(l.matrix)), biasOrder(l.biasOrder) {
        l.biasOrder = nullptr;
    }

    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    LUDecomposition<T, type, maxRow, maxColumn>& LUDecomposition<T, type, maxRow, maxColumn>::operator= (
            const LUDecomposition& l) {
        if(this != &l) {
            matrix = l.matrix;
            realloc(biasOrder, matrix.getRow() * sizeof(size_t));
        }
        return *this;
    }

    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    LUDecomposition<T, type, maxRow, maxColumn>& LUDecomposition<T, type, maxRow, maxColumn>::operator= (
            LUDecomposition&& l) noexcept {
        this->~LUDecomposition();
        matrix = std::move(l.matrix);
        biasOrder = l.biasOrder;
        l.biasOrder = nullptr;
    }

    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    LUDecomposition<T, type, maxRow, maxColumn>::~LUDecomposition() {
        free(biasOrder);
    }
    /*!
     * Apply LU Decomposition on a column of Matrix \from, save the result to Matrix \to.
     *
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009.32
     */
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    void LUDecomposition<T, type, maxRow, maxColumn>::decompositionColumn(size_t column) {
        const auto startAlphaIndex = column + 1;
        for (size_t j = 1; j < startAlphaIndex; ++j) {
            T temp(matrix(j, column));
            for (size_t k = 0; k < j; ++k)
                temp -= matrix(j, k) * matrix(k, column);
            matrix(j, column) = std::move(temp);
        }

        const auto r = matrix.getRow();
        for (size_t j = startAlphaIndex; j < r; ++j) {
            T temp(matrix(j, column));
            for (size_t k = 0; k < column; ++k)
                temp -= matrix(j, k) * matrix(k, column);
            matrix(j, column) = temp / matrix(column, column);
        }
    }
}

#endif