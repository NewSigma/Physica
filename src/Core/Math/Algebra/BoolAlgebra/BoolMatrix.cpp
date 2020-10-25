/*
 * Copyright 2020 WeiBo He.
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
#include "Physica/Core/Math/Algebra/BoolAlgebra/BoolMatrix.h"
#include <qglobal.h>

namespace Physica::Core {
    BoolMatrix::BoolMatrix(size_t column, size_t row) : arr(row) {
        Q_ASSERT(column > 0 && row > 0);
        for(size_t s = 0; s < row; ++s)
            arr.allocate(BitArray(column), s);
    }
    /*!
     * Return the multiple of two @class BoolMatrix: A(m * k) and B(k * n).
     *
     * Complexity: O(m * k * n)  possible to optimize.
     */
    BoolMatrix BoolMatrix::operator*(const BoolMatrix& m) const {
        const size_t c1 = getColumn();
        const size_t r1 = getRow();
        const size_t c2 = m.getColumn();
        Q_ASSERT(c1 == m.getRow());

        CStyleArray<BitArray, Dynamic> array(r1);
        for(size_t i = 0; i < r1; ++i) {
            BitArray row_i(c2);
            for(size_t j = 0; j < c2; ++j) {
                bool bit = false;
                for(size_t k = 0; k < c1; ++k) {
                    if((*this)[i][k] && m[k][j]) {
                        bit = true;
                        break;
                    }
                }
                row_i.setBit(bit, j);
            }
            array.allocate(std::move(row_i), i);
        }
        return BoolMatrix(std::move(array));
    }
    /*!
     * Complexity: O(row)
     */
    BoolMatrix BoolMatrix::operator&(const BoolMatrix& m) const {
        Q_ASSERT(hasSameSize(m));
        const size_t row = getRow();
        CStyleArray<BitArray, Dynamic> array(row);
        for(size_t s = 0; s < row; ++s)
            array.allocate(arr[s] & m.arr[s], s);
        return BoolMatrix(std::move(array));
    }
    /*!
     * Complexity: O(row)
     */
    BoolMatrix BoolMatrix::operator|(const BoolMatrix& m) const {
        Q_ASSERT(hasSameSize(m));
        const size_t row = getRow();
        CStyleArray<BitArray, Dynamic> array(row);
        for(size_t s = 0; s < row; ++s)
            array.allocate(arr[s] | m.arr[s], s);
        return BoolMatrix(std::move(array));
    }
    /*!
     * Complexity: O(row)
     */
    BoolMatrix BoolMatrix::operator~() const {
        const size_t row = getRow();
        CStyleArray<BitArray, Dynamic> array(row);
        for(size_t s = 0; s < row; ++s)
            array.allocate(~(arr[s]), s);
        return BoolMatrix(std::move(array));
    }
}