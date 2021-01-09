/*
 * Copyright 2020 WeiBo He.
 *
 * This file is part of Physica.

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
#ifndef PHYSICA_MATRIXOPERATION_H
#define PHYSICA_MATRIXOPERATION_H

namespace Physica::Core {
    /*!
     * This class contain static methods such as pivoting and eliminate.
     */
    template<class T, DenseMatrixType type, size_t maxRow, size_t maxColumn>
    class MatrixOperation {
    public:
        static inline size_t completePivoting(DenseMatrix<T, type, maxRow, maxColumn>& matrix, size_t column);
        static inline size_t partialPivoting(DenseMatrix<T, type, maxRow, maxColumn>& matrix, size_t column);
        static inline void impactPivoting(DenseMatrix<T, type, maxRow, maxColumn>& matrix, size_t row);
        static inline void upperEliminate(DenseMatrix<T, type, maxRow, maxColumn>& matrix, size_t index);
        static inline void lowerEliminate(DenseMatrix<T, type, maxRow, maxColumn>& matrix, size_t index);
        static inline void impactPivoting(DenseMatrix<T, type, maxRow, maxColumn>& matrix);
    };
    /*!
    * Select the main element of a column of Matrix and execute a column swap as well as a column swap to
    * make it stands at (column, column), return the origin column index of the main element.
    * The return values should be stored to recover the correct order of the solution.
    *
    * Complexity: O((rank - column) ^ 2)
    *
    * Reference:
    * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
    * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009.35
    */
    template<class T, DenseMatrixType type, size_t maxRow, size_t maxColumn>
    inline size_t MatrixOperation<T, type, maxRow, maxColumn>::completePivoting(
            DenseMatrix<T, type, maxRow, maxColumn>& matrix, size_t column) {
        const auto rank = matrix.getRow();
        Q_ASSERT(column < rank);
        size_t main_row_index = 0, main_column_index = 0;
        const T zero = T::getZero();
        const T* main = &zero;
        for(size_t i = column; i < rank; ++i) {
            for(size_t j = column; j < rank; ++j) {
                const auto* temp = &matrix(i, j);
                bool larger = absCompare(*main, *temp);
                main = larger ? main : temp;
                main_row_index = larger ? main_row_index : j;
                main_column_index = larger ? main_column_index : i;
            }
        }
        matrix.rowSwap(column, main_row_index);
        return main_column_index;
    }
    /*!
     * Select the main element of a column of Matrix and execute row swaps to make its row index equals to column index,
     * return the origin row index of the main element.
     *
     * Complexity: O(rank - column)
     *
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009.35
     */
    template<class T, DenseMatrixType type, size_t maxRow, size_t maxColumn>
    inline size_t MatrixOperation<T, type, maxRow, maxColumn>::partialPivoting(
            DenseMatrix<T, type, maxRow, maxColumn>& matrix, size_t column) {
        const auto rank = matrix.getRow();
        Q_ASSERT(column < rank);
        size_t main_index = column;
        const T* main = &matrix(column, column);
        for(size_t j = column + 1; j < rank; ++j) {
            const auto* temp = &matrix(j, column);
            bool larger = absCompare(*main, *temp);
            main = larger ? main : temp;
            main_index = larger ? main_index : j;
        }
        matrix.rowSwap(column, main_index);
        return main_index;
    }
    /*!
     * Divide the row by the element with the largest abstract value. \row is a row of a matrix.
     * \rank is the rank of the equations.
     *
     * Complexity: O(rank)
     *
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009.35
     */
    template<class T, DenseMatrixType type, size_t maxRow, size_t maxColumn>
    inline void MatrixOperation<T, type, maxRow, maxColumn>::impactPivoting(
            DenseMatrix<T, type, maxRow, maxColumn>& matrix, size_t row) {
        Q_ASSERT(row < matrix.getRow());
        const auto col = matrix.getColumn();
        const T* main = &matrix(row, 0);
        for(size_t i = 1; i < col; ++i) {
            const T* temp = &matrix(row, i);
            main = absCompare(*main, *temp) ? main : temp;
        }
        const T n = reciprocal(*main);
        for(size_t i = 0; i < col; ++i)
            matrix(row, i) *= n;
    }
    /*!
     * Eliminate elements at column index using row reduce, which are above the row index.
     * Warning: Element at (index, index) must not be zero. Or a divide by zero exception will be thrown.
     */
    template<class T, DenseMatrixType type, size_t maxRow, size_t maxColumn>
    inline void MatrixOperation<T, type, maxRow, maxColumn>::upperEliminate(
            DenseMatrix<T, type, maxRow, maxColumn>& matrix, size_t index) {
        for(size_t i = 0; i < index; ++i)
            matrix.rowReduce(index, i, index);
    }
    /*!
     * Eliminate elements at column index using row reduce, which are above the row index.
     * Warning: Element at (index, index) must not be zero. Or a divide by zero exception will be thrown.
     */
    template<class T, DenseMatrixType type, size_t maxRow, size_t maxColumn>
    inline void MatrixOperation<T, type, maxRow, maxColumn>::lowerEliminate(
            DenseMatrix<T, type, maxRow, maxColumn>& matrix, size_t index) {
        const auto r = matrix.getRow();
        for(size_t i = index + 1; i < r; ++i)
            matrix.rowReduce(index, i, index);
    }

    template<class T, DenseMatrixType type, size_t maxRow, size_t maxColumn>
    inline void MatrixOperation<T, type, maxRow, maxColumn>::impactPivoting(
            DenseMatrix<T, type, maxRow, maxColumn>& matrix) {
        const auto r = matrix.getRow();
        for(size_t i = 0; i < r; ++i)
            matrix.impactPivoting(i);
    }
}

#endif