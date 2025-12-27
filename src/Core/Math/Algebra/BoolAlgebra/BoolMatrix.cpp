/*
 * Copyright 2020-2025 Weibo He.
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

using namespace Physica;

BoolMatrix::BoolMatrix(size_t col, size_t row)
        : arr(row, BitArray(col)) {
    assert(col > 0 && row > 0);
}

BoolMatrix::BoolMatrix(size_t col, size_t row, bool initial)
        : arr(row, BitArray(col, initial)) {
    assert(col > 0 && row > 0);
}
/*!
 * Return the multiple of two @class BoolMatrix: A(m * k) and B(k * n).
 *
 * Complexity: O(m * k * n)  possible to optimize.
 */
BoolMatrix BoolMatrix::operator*(const BoolMatrix& m) const {
    assert(getCol() == m.getRow());
    return BoolMatrix(Array<BitArray>::generate(getRow(), [this, &m](size_t i) {
        const size_t c1 = getCol();
        const size_t c2 = m.getCol();
        BitArray row_i(c2);
        for (size_t j = 0; j < c2; ++j) {
            bool bit = false;
            for (size_t k = 0; k < c1; ++k) {
                if ((*this)[i][k] && m[k][j]) {
                    bit = true;
                    break;
                }
            }
            row_i.setBit(j, bit);
        }
        return row_i;
    }));
}
/*!
 * Complexity: O(row)
 */
BoolMatrix BoolMatrix::operator&(const BoolMatrix& m) const {
    assert(hasSameSize(m));
    return BoolMatrix(Array<BitArray>::generate(getRow(), [this, &m](size_t i) {
        return arr[i] & m.arr[i];
    }));
}
/*!
 * Complexity: O(row)
 */
BoolMatrix BoolMatrix::operator|(const BoolMatrix& m) const {
    assert(hasSameSize(m));
    return BoolMatrix(Array<BitArray>::generate(getRow(), [this, &m](size_t i) {
        return arr[i] | m.arr[i];
    }));
}
/*!
 * Complexity: O(row)
 */
BoolMatrix BoolMatrix::operator~() const {
    return BoolMatrix(Array<BitArray>::generate(getRow(), [this](size_t i) {
        return ~(arr[i]);
    }));
}
