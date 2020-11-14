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
#ifndef PHYSICA_BOOLMATRIX_H
#define PHYSICA_BOOLMATRIX_H

#include "Physica/Utils/BitArray.h"
#include "Physica/Utils/Container/CStyleArray/CStyleArray.h"

namespace Physica::Core {
    using namespace Utils;
    /*!
     * Optimize:
     * 1. Fixed size matrix
     * 2. Implement a column matrix.
     */
    class BoolMatrix {
        CStyleArray<BitArray, Dynamic> arr;
    public:
        BoolMatrix(size_t column, size_t row);
        BoolMatrix(const BoolMatrix& m) = default;
        BoolMatrix(BoolMatrix&& m) noexcept : arr(std::move(m.arr)) {}
        ~BoolMatrix() = default;
        /* Operators */
        BoolMatrix& operator=(const BoolMatrix& m) { if(this != &m) arr = m.arr; return *this; }
        BoolMatrix& operator=(BoolMatrix&& m) noexcept { arr = std::move(m.arr); return *this; }
        BitArray& operator[](size_t s) { return arr[s]; }
        const BitArray& operator[](size_t s) const { return arr[s]; }
        [[nodiscard]] BoolMatrix operator*(const BoolMatrix& m) const;
        [[nodiscard]] BoolMatrix operator&(const BoolMatrix& m) const;
        [[nodiscard]] BoolMatrix operator|(const BoolMatrix& m) const;
        [[nodiscard]] BoolMatrix operator~() const;
        /* Operations */
        [[nodiscard]] inline bool hasSameSize(const BoolMatrix& m) const;
        /* Getters */
        [[nodiscard]] size_t getColumn() const { return arr[0].getLength(); }
        [[nodiscard]] size_t getRow() const { return arr.getLength(); }
    private:
        /*!
         * Construct a BoolMatrix from its members, declared private to avoid improper uses.
         */
        explicit BoolMatrix(CStyleArray<BitArray, Dynamic>&& arr) : arr(std::move(arr)) {}
    };

    inline bool BoolMatrix::hasSameSize(const BoolMatrix& m) const {
        return getRow() == m.getRow() && getColumn() == m.getColumn();
    }
}

#endif
