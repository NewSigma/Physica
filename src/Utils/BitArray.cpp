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
#include "Physica/Utils/BitArray.h"

namespace Physica::Utils {
    /*!
     * bitCount is the number of bool types you need.
     * bitCount is different from BitArray::length.
     */
    BitArray::BitArray(size_t bitCount) : bitCount(bitCount), arr(new unsigned char[getLength()]) {}

    BitArray::BitArray(const BitArray& array)
            : bitCount(array.bitCount)
            , arr(new unsigned char[array.getLength()]) {
        const size_t length = getLength();
        for(unsigned int i = 0; i < length; ++i)
            arr[i] = array.arr[i];
    }

    BitArray::BitArray(BitArray&& array) noexcept : arr(array.arr), bitCount(array.bitCount) {
        array.arr = nullptr;
    }

    BitArray::~BitArray() {
        delete[] arr;
    }

    BitArray& BitArray::operator=(const BitArray& array) {
        if(this != &array) {
            this->~BitArray();
            bitCount = array.bitCount;
            const size_t length = getLength();
            arr = reinterpret_cast<unsigned char*>(malloc(length * sizeof(unsigned char)));
            for(unsigned int i = 0; i < length; ++i)
                arr[i] = array.arr[i];
        }
        return *this;
    }

    BitArray& BitArray::operator=(BitArray&& array) noexcept {
        if(this != &array)
            this->~BitArray();
        bitCount = array.bitCount;
        arr = array.arr;
        return *this;
    }
    /*!
     * Access the s.th bool, s starts from 0.
     */
    bool BitArray::operator[](size_t s) const {
        Q_ASSERT(bitCount > s);
        s += 1;
        const auto quotient = s / CHAR_BIT;
        const auto reminder = s - quotient * CHAR_BIT;
        if(reminder)
            return arr[quotient] & (1U << (reminder - 1U));
        return arr[quotient - 1] & (1U << (CHAR_BIT - 1U));
    }

    BitArray BitArray::operator&(const BitArray& array) const {
        Q_ASSERT(bitCount == array.bitCount);
        const size_t length = getLength();
        auto newArr = reinterpret_cast<unsigned char*>(malloc(length * sizeof(unsigned char)));
        for(size_t s = 0; s < length; ++s)
            newArr[s] = arr[s] & array.arr[s];
        return BitArray(newArr, length);
    }

    BitArray BitArray::operator|(const BitArray& array) const {
        Q_ASSERT(bitCount == array.bitCount);
        const size_t length = getLength();
        auto newArr = reinterpret_cast<unsigned char*>(malloc(length * sizeof(unsigned char)));
        for(size_t s = 0; s < length; ++s)
            newArr[s] = arr[s] | array.arr[s];
        return BitArray(newArr, length);
    }

    BitArray BitArray::operator~() const {
        const size_t length = getLength();
        auto newArr = reinterpret_cast<unsigned char*>(malloc(length * sizeof(unsigned char)));
        for(size_t s = 0; s < length; ++s)
            newArr[s] = ~arr[s];
        return BitArray(newArr, length);
    }

    void BitArray::setBit(bool b, size_t s) {
        Q_ASSERT(bitCount > s);
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
    BitArray::BitArray(unsigned char* arr, size_t bitCount) : arr(arr), bitCount(bitCount) {}
}