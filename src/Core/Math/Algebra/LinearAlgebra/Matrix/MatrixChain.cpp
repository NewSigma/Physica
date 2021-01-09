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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixChain.h"

namespace Physica::Core {
    /*!
     * \class MatrixChain provides a effective method to calculate the production of a chain of matrices.
     *
     * Note: matrices should be initialized and pass into the MatrixChain before calling solve(),
     * matrices must not be deleted before the calculation finished.
     *
     * Warning:
     * 1.If the chain is too long, size_t or stack may overflow and lead to wrong results.
     * 2.length should be at least 2 or solve() will simply return the matrix passed in
     * , which may result in double free.
     *
     * Reference: 算法导论 第三版 Page: 210-215
     */
    template<class T, DenseMatrixType type, size_t maxRow, size_t maxColumn>
    MatrixChain<T, type, maxRow, maxColumn>::MatrixChain(size_t length)
            : chain(new DenseMatrix<T, type, maxRow, maxColumn>*[length])
            , length(length), price(new size_t*[length]), point(new size_t*[length - 1]) {
        const auto length_1 = length - 1;
        for(size_t i = 0; i < length_1; ++i) {
            price[i] = new size_t[length];
            point[i] = new size_t[length_1];
        }
        //price is longer than point by 1.
        price[length_1] = new size_t[length];
    }

    template<class T, DenseMatrixType type, size_t maxRow, size_t maxColumn>
    MatrixChain<T, type, maxRow, maxColumn>::~MatrixChain() {
        Q_UNUSED(type)
        Q_UNUSED(maxRow)
        Q_UNUSED(maxColumn)
        delete[] chain;
        const auto length_1 = length - 1;
        for(size_t i = 0; i < length_1; ++i) {
            delete[] price[i];
            delete[] point[i];
        }
        //price is longer than point by 1.
        delete[] price[length_1];
        delete[] price;
        delete[] point;
    }
    //!Optimize: Only half of the space of price and point is used. Maybe change them into a 1D array.
    template<class T, DenseMatrixType type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, Dynamic, Dynamic> MatrixChain<T, type, maxRow, maxColumn>::solve() {
        Q_UNUSED(type)
        Q_UNUSED(maxRow)
        Q_UNUSED(maxColumn)
        for(size_t i = 0; i < length; ++i)
            price[i][i] = 0;

        for(size_t l = 1; l < length; ++l) { //(l + 1) is the length of sub-chain.
            const auto m_max = length - l;
            for(size_t m = 0; m < m_max; ++m) { //m is the start index of each sub-chain.
                const auto m_end = m + l; //End index of each sub-chain.
                auto n = m; //Cut (m ... ) into (m ... n) and (n + 1 ... m_end).
                auto chain_n = chain[n];
                /* Handle n = m */ {
                    price[m][m_end] = price[n + 1][m_end] + chain_n->getRow() * chain_n->getColumn() * chain[n + 1]->getRow();
                    point[m][m_end] = n;
                    ++n;
                }
                for(; n < m_end; ++n) {
                    chain_n = chain[n];
                    size_t temp = price[m][n] + price[n + 1][m_end]
                                                + chain_n->getRow() * chain_n->getColumn() * chain[n + 1]->getRow();
                    price[m][m_end] = temp < price[m][m_end] ? temp : price[m][m_end];
                    point[m][m_end] = temp < price[m][m_end] ? n : point[m][m_end];
                }
                /* Handle n = m_end */ {
                    chain_n = chain[n];
                    size_t temp = price[m][n] + chain_n->getRow() * chain_n->getColumn() * chain[n + 1]->getRow();
                    price[m][m_end] = temp < price[m][m_end] ? temp : price[m][m_end];
                    point[m][m_end] = temp < price[m][m_end] ? n : point[m][m_end];
                }
            }
        }
        return multiply(0, length - 1);
    }
    //!Both \from and \to are included.
    template<class T, DenseMatrixType type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, Dynamic, Dynamic> MatrixChain<T, type, maxRow, maxColumn>::multiply(size_t from, size_t to) {
        Q_UNUSED(maxRow)
        Q_UNUSED(maxColumn)
        if(from == to)
            return DenseMatrix<T, type, Dynamic, Dynamic>(*chain[from]);
        const auto cut_at = point[from][to];
        auto first = multiply(from, cut_at);
        auto second = multiply(cut_at + 1, to);
        return first * second;
    }
}