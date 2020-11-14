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
#ifndef PHYSICA_BITARRAY_H
#define PHYSICA_BITARRAY_H

#include <cstddef>

namespace Physica::Utils {
    /*!
     * BitArray is similar to bool[] but it has a better use of space.
     */
    class BitArray {
        //Every bit of arr stands for a bool.
        unsigned char* __restrict__ arr;
        size_t bitCount;
    public:
        explicit BitArray(size_t bitCount);
        BitArray(const BitArray& array);
        BitArray(BitArray&& array) noexcept;
        ~BitArray();
        /* Operators */
        BitArray& operator=(const BitArray& array);
        BitArray& operator=(BitArray&& array) noexcept;
        bool operator[](size_t s) const;
        BitArray operator&(const BitArray& array) const;
        BitArray operator|(const BitArray& array) const;
        BitArray operator~() const;
        /* Operations */
        void setBit(bool b, size_t s);
        /* Getters */
        [[nodiscard]] inline size_t getLength() const;
    private:
        BitArray(unsigned char* arr, size_t bitCount);
    };

    inline size_t BitArray::getLength() const {
        //((bitCount >> 3U) + (bitCount & 8U)) is the upper approximation of (bitCount / 8).
        return (bitCount >> 3U) + (bitCount & 8U);
    }
}

#endif
