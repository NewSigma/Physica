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
#include <cstdlib>
#include <climits>
#include "qglobal.h"
#include "Physica/Core/Utils/BitArray.h"

namespace Physica::Core {
    /*!
     * bitCount is the number of bool types you need.
     * bitCount is different from BitArray::length.
     */
    BitArray::BitArray(unsigned int bitCount) : length((bitCount >> 3U) + (bitCount & 8U)) { //((bitCount >> 3U) + (bitCount & 8U)) is the upper approximation of (bitCount / 8).
        arr = reinterpret_cast<unsigned char*>(malloc(length * sizeof(unsigned char)));
    }

    BitArray::BitArray(const BitArray& array)
            : arr(arr = reinterpret_cast<unsigned char*>(malloc(array.length * sizeof(unsigned char)))), length(array.length) {
        for(unsigned int i = 0; i < length; ++i)
            arr[i] = array.arr[i];
    }

    BitArray::BitArray(BitArray&& array) noexcept : arr(array.arr), length(array.length) {
        if(this != &array)
            array.arr = nullptr;
    }

    BitArray::~BitArray() {
        delete[] arr;
    }

    BitArray& BitArray::operator=(const BitArray& array) {
        if(this != &array) {
            this->~BitArray();
            length = array.length;
            arr = reinterpret_cast<unsigned char*>(malloc(length * sizeof(unsigned char)));
            for(unsigned int i = 0; i < length; ++i)
                arr[i] = array.arr[i];
        }
        return *this;
    }

    BitArray& BitArray::operator=(BitArray&& array) noexcept {
        if(this != &array)
            this->~BitArray();
        length = array.length;
        arr = array.arr;
        return *this;
    }
    /*!
     * Access the s.th bool, s starts from 0.
     */
    bool BitArray::operator[](size_t s) const {
        Q_ASSERT(length * CHAR_BIT > s);
        s += 1;
        const auto quotient = s / CHAR_BIT;
        const auto reminder = s - quotient * CHAR_BIT;
        if(reminder)
            return arr[quotient] & (1U << (reminder - 1U));
        return arr[quotient - 1] & (1U << (CHAR_BIT - 1U));
    }

    BitArray BitArray::operator&(const BitArray& array) const {
        auto minLength = std::min(length, array.length);
        auto newArr = reinterpret_cast<unsigned char*>(malloc(minLength * sizeof(unsigned char)));
        for(size_t s = 0; s < minLength; ++s)
            newArr[s] = arr[s] & array.arr[s];
        return BitArray(newArr, minLength);
    }

    BitArray BitArray::operator|(const BitArray& array) const {
        auto minLength = std::min(length, array.length);
        auto newArr = reinterpret_cast<unsigned char*>(malloc(minLength * sizeof(unsigned char)));
        for(size_t s = 0; s < minLength; ++s)
            newArr[s] = arr[s] | array.arr[s];
        return BitArray(newArr, minLength);
    }

    BitArray BitArray::operator~() const {
        auto newArr = reinterpret_cast<unsigned char*>(malloc(length * sizeof(unsigned char)));
        for(size_t s = 0; s < length; ++s)
            newArr[s] = ~arr[s];
        return BitArray(newArr, length);
    }

    void BitArray::setBit(bool b, size_t s) {
        Q_ASSERT(length * CHAR_BIT > s);
        s += 1;
        const auto quotient = s / CHAR_BIT;
        const auto reminder = s - quotient * CHAR_BIT;
        if(b) {
            if(reminder)
                arr[quotient] |= (1U << (reminder - 1U));
            else
                arr[quotient - 1] |= (1U << (CHAR_BIT - 1U));
        }
        else {
            if(reminder)
                arr[quotient] &= ~(1U << (reminder - 1U));
            else
                arr[quotient - 1] &= ~(1U << (CHAR_BIT - 1U));
        }
    }
    /*!
     * Construct a BitArray directly from its members. Declared private to avoid incorrect uses.
     * Length of arr must be @param length.
     */
    BitArray::BitArray(unsigned char* arr, size_t length) : arr(arr), length(length) {}
}