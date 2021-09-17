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
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "Physica/Core/MultiPrecision/Integer.h"
#include "Physica/Core/MultiPrecision/BasicImpl/Util/Bitwise.h"
#include "Physica/Core/MultiPrecision/BasicImpl/Convert.h"

namespace Physica::Core {
    Integer::Integer(int i)
            : byte(reinterpret_cast<MPUnit*>(malloc(sizeof(MPUnit))))
            , length(i >= 0 ? 1 : -1) {
        byte[0] = i >= 0 ? i : -i;
    }

    Integer::Integer(const Integer& toCopy)
            : byte(reinterpret_cast<MPUnit*>(malloc(toCopy.getSize() * sizeof(MPUnit))))
            , length(toCopy.length) {
        memcpy(byte, toCopy.byte, getSize() * sizeof(MPUnit));
    }

    Integer::Integer(Integer&& toMove) noexcept
            : byte(toMove.byte)
            , length(toMove.length) {
        toMove.byte = nullptr;
    }

    Integer::~Integer() {
        free(byte);
    }

    Integer& Integer::operator=(const Integer& toCopy) {
        if (this != &toCopy) {
            this->~Integer();
            length = toCopy.length;
            int size = getSize();
            byte = reinterpret_cast<MPUnit*>(malloc(size * sizeof(MPUnit)));
            memcpy(byte, toCopy.byte, size * sizeof(MPUnit));
        }
        return *this;
    }

    Integer& Integer::operator= (Integer&& toMove) noexcept {
        this->~Integer();
        byte = toMove.byte;
        toMove.byte = nullptr;
        length = toMove.length;
        return *this;
    }

    Integer Integer::operator+ (const Integer& i) const {
        if (isZero())
            return i;
        if (i.isZero())
            return *this;
        if (!matchSign(*this, i))
            return length > 0 ? Integer::integerSubImpl(*this, i) : integerSubImpl(i, *this);
        Integer abs_i1(byte, length > 0 ? length : -length);
        Integer abs_i2(i.byte, i.length > 0 ? i.length : -i.length);
        Integer result = integerAddImpl(abs_i1, abs_i2);
        abs_i1.byte = nullptr;
        abs_i2.byte = nullptr;
        return length > 0 ? result : -result;
    }

    Integer Integer::operator-(const Integer& i) const {
        if (isZero())
            return -i;
        if (i.isZero())
            return -*this;
        if (length > 0) {
            if (i.length > 0)
                return integerSubImpl(*this, i);
            else {
                Integer abs_i2(i.byte, -i.length);
                Integer result = integerAddImpl(*this, abs_i2);
                abs_i2.byte = nullptr;
                return result;
            }
        }
        else {
            if (i.length > 0) {
                Integer abs_i1(byte, -length);
                Integer result = integerAddImpl(abs_i1, i);
                abs_i1.byte = nullptr;
                return -result;
            }
            else {
                Integer abs_i1(byte, -length);
                Integer abs_i2(i.byte, -i.length);
                Integer result = integerSubImpl(abs_i2, abs_i1);
                abs_i1.byte = nullptr;
                abs_i2.byte = nullptr;
                return result;
            }
        }
    }

    Integer Integer::operator*(const Integer& i) const {
        if (isZero() || i.isZero())
            return 0;
        if (*this == 1)
            return i;
        if (i == 1)
            return *this;
        const int size1 = getSize();
        const int size2 = i.getSize();
        //Estimate the ed of result first. we will calculate it accurately later.
        auto resultLength = size1 + size2;
        auto* __restrict resultByte = reinterpret_cast<MPUnit*>(calloc(resultLength, sizeof(MPUnit)));
        for (int j = 0; j < size1; ++j)
            resultByte[j + size2] = mulAddArrByWord(resultByte + j, i.byte, size2, byte[j]);
        if (resultByte[resultLength - 1] == 0) {
            --resultLength;
            resultByte = reinterpret_cast<MPUnit*>(realloc(resultByte, resultLength * sizeof(MPUnit)));
        }
        return Integer(resultByte, matchSign(*this, i) ? resultLength : -resultLength);
    }

    Integer Integer::operator/(const Integer& i) const {
        if (i.isZero())
            throw DivideByZeroException();
        if (isZero())
            return 0;
        if (i == 1)
            return *this;
        const int i1_size = getSize();
        const int i2_size = i.getSize();

        auto arr1_len = std::max(i1_size, i2_size) + 1;
        auto i1_blank = arr1_len - i1_size;
        auto arr1 = new MPUnit[arr1_len];
        memcpy(arr1 + i1_blank, byte, i1_size * sizeof(MPUnit));
        memset(arr1, 0, i1_blank * sizeof(MPUnit));
        //Size of arr2 is arranged 1 less than arr1.
        auto arr2_len = arr1_len - 1;
        auto i2_blank = arr2_len - i2_size;
        auto arr2 = new MPUnit[arr2_len];
        memcpy(arr2 + i2_blank, i.byte, i2_size * sizeof(MPUnit));
        memset(arr2, 0, i2_blank * sizeof(MPUnit));
        /*
         * We shift s1 and s2, making the less highest bit of s1 is set and the highest bit of s2 is set
         * to meet the acquirement of the function divArrByFullArrWith1Word().
         */
        const int i1_shift = static_cast<int>(countLeadingZeros(byte[i1_size - 1])) - 1;
        if(i1_shift > 0)
            byteLeftShiftEq(arr1, arr1_len, i1_shift);
        else
            byteRightShiftEq(arr1, arr1_len, -i1_shift);
        const int i2_shift = static_cast<int>(countLeadingZeros(i.byte[i2_size - 1]));
        byteLeftShiftEq(arr2, arr2_len, i2_shift);
        ////////////////////////////////Calculate cursory first//////////////////////////////////////
        //Estimate the length of result.
        int resultLength = arr2_len;
        auto* __restrict resultByte = reinterpret_cast<MPUnit*>(malloc(resultLength * sizeof(MPUnit)));
        for(int j = resultLength - 1; j >= 0; --j) {
            resultByte[j] = divArrByFullArrWith1Word(arr1, arr2, arr2_len);
            arr1[arr2_len] -= mulSubArrByWord(arr1, arr2, arr2_len, byte[j]);
            byteLeftShiftEq(arr1, arr1_len, MPUnitWidth);
        }
        delete[] arr1;
        delete[] arr2;
        ////////////////////////////////////Out put////////////////////////////////////////
        Integer temp(resultByte, matchSign(*this, i) ? resultLength : -resultLength);
        return temp >> ((i1_blank - i2_blank) * static_cast<int>(MPUnitWidth) + (i1_shift - i2_shift));
    }

    Integer Integer::operator%(const Integer& i) const {
        return *this - *this / i * i;
    }

    Integer Integer::operator<<(int bits) const {
        if (bits == 0)
            return *this;
        if (bits < 0)
            return *this >> -bits;
        const int size = getSize();
        const int quotient = bits / MPUnitWidth; //NOLINT: quotient < INT_MAX
        const unsigned int remainder = bits - quotient * MPUnitWidth;
        //If remainder = 0, we must return directly because shifting a MPUnit for MPUnitWidth bits is a undefined behavior.
        MPUnit* __restrict resultByte;
        int resultLength = size + quotient;
        if (remainder == 0) {
            resultByte = reinterpret_cast<MPUnit*>(malloc(resultLength * sizeof(MPUnit)));
            memset(resultByte, 0, quotient * sizeof(MPUnit));
            memcpy(resultByte + quotient, byte, resultLength * sizeof(MPUnit));
            return Integer(resultByte, resultLength);
        }
        const bool carry = countLeadingZeros(byte[size - 1]) < remainder;
        resultLength += carry;
        resultByte = reinterpret_cast<MPUnit*>(malloc(resultLength * sizeof(MPUnit)));
        memset(resultByte, 0, quotient * sizeof(MPUnit));
        resultByte[quotient] = 0;
        const int size_1 = quotient + size - 1;
        for(int i = quotient; i < size_1; ++i) {
            resultByte[i] |= byte[i] << remainder;
            resultByte[i + 1] = byte[i] >> (MPUnitWidth - remainder);
        }
        resultByte[size_1] |= byte[size_1] << remainder;
        if(carry)
            resultByte[size + quotient] = byte[size_1] >> (MPUnitWidth - remainder);
        return Integer(resultByte, resultLength);
    }

    Integer Integer::operator>>(int bits) const {
        if(bits == 0)
            return *this;
        if(bits < 0)
            return *this << -bits;
        const int size = getSize();
        const int quotient = bits / MPUnitWidth; //NOLINT: quotient < INT_MAX
        const unsigned int remainder = bits - quotient * MPUnitWidth;
        //If remainder = 0, we must return directly because shifting a MPUnit for MPUnitWidth bits is a undefined behavior.
        MPUnit* __restrict resultByte;
        int resultLength = size > quotient ? size - quotient : 1;
        if (remainder == 0) {
            resultByte = reinterpret_cast<MPUnit*>(malloc(resultLength * sizeof(MPUnit)));
            memcpy(resultByte, byte + quotient, resultLength * sizeof(MPUnit));
            resultByte[0] = size > quotient ? resultByte[0] : 0;
            return Integer(resultByte, resultLength);
        }
        const bool carry = (countLeadingZeros(byte[size - 1]) + remainder) >= MPUnitWidth;
        resultLength -= carry;
        resultByte = reinterpret_cast<MPUnit*>(malloc(resultLength * sizeof(MPUnit)));
        if(!carry)
            resultByte[size - 1] = byte[size - 1] >> remainder;

        for(int i = size - 2; i >= 0; --i) {
            resultByte[i] = byte[i] >> remainder;
            resultByte[i] |= byte[i + 1] << (MPUnitWidth - remainder);
        }
        Integer result(resultByte, resultLength);
        result.cutZero();
        return result;
    }

    Integer Integer::operator-() const {
        Integer copy(*this);
        copy.length = -copy.length;
        return copy;
    }

    bool Integer::operator>(const Integer& i) const {
        if (length == i.length) {
            for (int j = getSize() - 1; j >= 0; --j)
                if (byte[j] > i.byte[j])
                    return true;
            return false;
        }
        return length > i.length;
    }

    bool Integer::operator<(const Integer& i) const {
        if (length == i.length) {
            for (int j = getSize() - 1; j >= 0; --j)
                if (byte[j] < i.byte[j])
                    return true;
            return false;
        }
        return length < i.length;
    }

    bool Integer::operator==(const Integer& i) const {
        if (length != i.length)
            return false;
        for (int j = 0; j < length; ++j)
            if (byte[j] != i.byte[j])
                return false;
        return true;
    }

    Integer& Integer::operator++() {
        *this += 1;
        return *this;
    }

    Integer& Integer::operator--() {
        *this -= 1;
        return *this;
    }

    Integer Integer::operator++(int) {
        Integer temp(*this);
        *this += 1;
        return temp;
    }

    Integer Integer::operator--(int) {
        Integer temp(*this);
        *this -= 1;
        return temp;
    }

    Integer::operator double() const {
        if (isZero())
            return 0;
        return Internal::convertDoubleImpl(length, length - 1, byte);
    }

    void Integer::swap(Integer& i) noexcept {
        std::swap(byte, i.byte);
        std::swap(length, i.length);
    }
    /**
     * return true if the abstract value of i1 is larger or equal than the abstract value of i2.
     * return false if the abstract value of i1 is smaller to the abstract value of i2.
     */
    bool Integer::absCompare(const Integer& i1, const Integer& i2) {
        if (i1.length == i2.length) {
            for (int i = i1.length - 1; i >= 0; --i) {
                if (i1.byte[i] == i2.byte[i])
                    continue;
                return i1.byte[i] > i2.byte[i];
            }
        }
        return std::abs(i1.length) > std::abs(i2.length);
    }

    void Integer::cutZero() {
        if (!isZero()) {
            for (int i = length - 1; i > 0; --i) {
                if (byte[i] != 0)
                    break;
                --length;
            }
            byte = reinterpret_cast<MPUnit*>(realloc(byte, length * sizeof(MPUnit)));
        }
    }
}