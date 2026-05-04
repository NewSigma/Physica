/*
 * Copyright 2026 Weibo He.
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
#include "Physica/Core/Scalar/ScalarImpl/ArrayArithmetic.h"
#include "Physica/Core/Utils/Allocator/HostAllocator.h"

using namespace Physica;

MPUnit Physica::addArrWithArr(MPUnit* __restrict result, const MPUnit* __restrict from, const MPUnit* __restrict to, size_t len) noexcept {
    MPUnit carry = 0;
    for (size_t i = 0; i < len; ++i) {
        MPUnit from_i = from[i];
        MPUnit temp = from_i + to[i] + carry;
        result[i] = temp;
        carry = temp < from_i;
    }
    return carry;
}

MPUnit Physica::addArrWithArrEq(const MPUnit* __restrict from, MPUnit* __restrict to, size_t len) noexcept {
    MPUnit carry = 0;
    for (size_t i = 0; i < len; ++i) {
        MPUnit from_i = from[i];
        MPUnit add = from_i + to[i] + carry;
        to[i] = add;
        carry = add < from_i;
    }
    return carry;
}

void Physica::sub2WordByWord(MPUnit& high, MPUnit& low, MPUnit n1_high, MPUnit n1_low, MPUnit n2) noexcept {
    low = n1_low - n2;
    high = n1_high - low > n1_low;
}

void Physica::sub2WordByWord(MPUnit& high, MPUnit& low, MPUnit n2) noexcept {
    // Use n2 as a temp variable, instead of creating a new temp variable.
    n2 = low - n2;
    high -= low > n2;
    low = n2;
}

bool Physica::subArrByArr(MPUnit* __restrict result, const MPUnit* __restrict arr1, const MPUnit* __restrict arr2, size_t len) noexcept {
    MPUnit carry = 0, temp1, temp2;
    bool pre_carry;
    for (size_t i = 0; i < len; ++i) {
        temp1 = arr1[i];
        temp2 = temp1 - arr2[i];
        pre_carry = temp1 < temp2;
        temp1 = temp2 - carry;
        result[i] = temp1;
        carry = MPUnit(pre_carry | (temp1 > temp2));
    }
    return bool(carry);
}

bool Physica::subArrByArrEq(MPUnit* __restrict arr1, const MPUnit* __restrict arr2, size_t len) noexcept {
    MPUnit carry = 0, temp1, temp2;
    bool pre_carry;
    for (size_t i = 0; i < len; ++i) {
        temp1 = arr1[i];
        temp2 = temp1 - arr2[i];
        pre_carry = temp1 < temp2;
        temp1 = temp2 - carry;
        arr1[i] = temp1;
        carry = MPUnit(pre_carry | (temp1 > temp2));
    }
    return bool(carry);
}

MPUnit Physica::mulWordByWordHigh(MPUnit n1, MPUnit n2) noexcept {
    if constexpr (!IsMSVC()) {
#ifdef __GNUC__
        MPUnit result;
        if constexpr (PhysicaWordSize == 64) {
            asm(
                    "mulq %2\n\t"
                    "movq %%rdx, %0"
                    : "=a"(result)
                    : "a"(n1), "rm"(n2)
                    : "%rdx");
        }
        else {
            asm(
                    "mull %2"
                    : "=d"(result)
                    : "a"(n1), "rm"(n2));
        }
        return result;
#else
        noImpl();
#endif
    }
    else {
        MPUnit n1_low = n1 & MPUnitLowerMask;
        MPUnit n1_high = n1 >> (64U / 2U);
        MPUnit n2_low = n2 & MPUnitLowerMask;
        MPUnit n2_high = n2 >> (64U / 2U);

        MPUnit ll = n1_low * n2_low;
        MPUnit lh = n1_low * n2_high;
        MPUnit hl = n1_high * n2_low;
        MPUnit hh = n1_high * n2_high;

        lh += ll >> (64U / 2U);
        lh += hl;
        hh += static_cast<MPUnit>(lh < hl) << (64U / 2U);
        return hh + (lh >> (64U / 2U));
    }
}

void Physica::mulWordByWord(MPUnit& high, MPUnit& low, MPUnit n1, MPUnit n2) noexcept {
    if constexpr (!IsMSVC()) {
#ifdef __GNUC__
        if constexpr (PhysicaWordSize == 64) {
            asm(
                    "mulq %3"
                    : "=d"(high), "=a"(low)
                    : "a"(n1), "rm"(n2));
        }
        else {
            asm(
                    "mull %3"
                    : "=d"(high), "=a"(low)
                    : "a"(n1), "rm"(n2));
        }
#else
        noImpl();
#endif
    }
    else {
        MPUnit n1_low = n1 & MPUnitLowerMask;
        MPUnit n1_high = n1 >> (MPUnitWidth / 2U);
        MPUnit n2_low = n2 & MPUnitLowerMask;
        MPUnit n2_high = n2 >> (MPUnitWidth / 2U);

        auto ll = n1_low * n2_low;
        auto lh = n1_low * n2_high;
        auto hl = n1_high * n2_low;
        auto hh = n1_high * n2_high;

        lh += ll >> (MPUnitWidth / 2U);
        lh += hl;
        hh += static_cast<MPUnit>(lh < hl) << (MPUnitWidth / 2U);
        high = hh + (lh >> (MPUnitWidth / 2U));
        low = (lh << (MPUnitWidth / 2U)) + (ll & MPUnitLowerMask);
    }
}

MPUnit Physica::mulArrByWord(MPUnit* __restrict result, const MPUnit* __restrict arr, size_t length, MPUnit n) noexcept {
    MPUnit carry = 0, high, low;
    for (size_t i = 0; i < length; ++i) {
        mulWordByWord(high, low, arr[i], n);
        low += carry;
        carry = (low < carry) + high;
        result[i] = low;
    }
    return carry;
}

MPUnit Physica::mulAddArrByWord(MPUnit* __restrict result, const MPUnit* __restrict arr, size_t length, MPUnit n) noexcept {
    MPUnit carry = 0, high, low;
    for (size_t i = 0; i < length; ++i) {
        mulWordByWord(high, low, arr[i], n);
        low += carry;
        carry = (low < carry) + high;
        result[i] += low;
        carry += result[i] < low;
    }
    return carry;
}

MPUnit Physica::mulSubArrByWord(MPUnit* __restrict result, const MPUnit* __restrict arr, size_t length, MPUnit n) noexcept {
    MPUnit carry = 0, high, low;
    for (size_t i = 0; i < length; ++i) {
        mulWordByWord(high, low, arr[i], n);
        low = result[i] - low;
        high += result[i] < low;
        result[i] = low;
        low -= carry;
        carry = high + (result[i] < low);
        result[i] = low;
    }
    return carry;
}

MPUnit Physica::getInverse(MPUnit unit) noexcept {
    MPUnit unit0 = unit & 1U;
    if constexpr (PhysicaWordSize == 64) {
        MPUnit unit9 = unit >> 55U;
        MPUnit unit40 = (unit >> 24U) + 1;
        MPUnit unit63 = (unit >> 1U) + unit0;
        MPUnit v0 = ((static_cast<MPUnit>(1) << 19U) - 3 * (static_cast<MPUnit>(1) << 8U)) / unit9;
        MPUnit v1 = (static_cast<MPUnit>(1) << 11U) * v0 - ((v0 * v0 * unit40) >> 40U) - 1;
        MPUnit v2 = (static_cast<MPUnit>(1) << 13U) * v1
                  + ((((static_cast<MPUnit>(1) << 60U) - v1 * unit40) * v1) >> 47U);
        MPUnit e = (v2 >> 1U) * unit0 - v2 * unit63;
        MPUnit v3 = (static_cast<MPUnit>(1) << 31U) * v2 + (mulWordByWordHigh(v2, e) >> 1U);
        MPUnit high, low;
        mulWordByWord(high, low, v3, unit);
        MPUnit v4 = v3 - unit - high - (low + unit < low);
        return v4;
    }
    else {
        MPUnit unit10 = unit >> 22U;
        MPUnit unit21 = (unit >> 11U) + 1;
        MPUnit unit31 = (unit >> 1U) + unit10;
        MPUnit v0 = ((static_cast<MPUnit>(1) << 24U) - (static_cast<MPUnit>(1) << 14U)
                     + (static_cast<MPUnit>(1) << 9U))
                  / unit10;
        MPUnit v1 = (static_cast<MPUnit>(1) << 4U) * v0 - mulWordByWordHigh(v0 * v0, unit21) - 1;
        MPUnit e = (v1 >> 1U) * unit0 - v1 * unit31;
        MPUnit v2 = (static_cast<MPUnit>(1) << 15U) * v1 - (mulWordByWordHigh(v1, e) >> 1U);
        MPUnit high, low;
        mulWordByWord(high, low, v2, unit);
        MPUnit v3 = v2 - unit - high - (low + unit < low);
        return v3;
    }
}

void Physica::div2WordByFullWord(MPUnit& quotient, MPUnit& remainder, MPUnit high, MPUnit low, MPUnit divisor) noexcept {
    assert(high < divisor);
    assert(divisor & MPUnitHighestBitMask);
    MPUnit quotient2;
    mulWordByWord(quotient, quotient2, high, getInverse(divisor));

    auto temp = quotient2;
    quotient2 += low;
    quotient += high + (temp > quotient2) + 1;
    remainder = low - quotient * divisor;
    if (remainder > quotient2) {
        --quotient;
        remainder += divisor;
    }
    if (remainder >= divisor) {
        ++quotient;
        remainder -= divisor;
    }
}

MPUnit Physica::div2WordByFullWordQ(MPUnit high, MPUnit low, MPUnit divisor) noexcept {
    assert(high < divisor);
    assert(divisor & MPUnitHighestBitMask);
    MPUnit quotient, quotient2;
    mulWordByWord(quotient, quotient2, high, getInverse(divisor));

    auto temp = quotient2;
    quotient2 += low;
    quotient += high + (temp > quotient2) + 1;
    MPUnit remainder = low - quotient * divisor;
    if (remainder > quotient2) {
        --quotient;
        remainder += divisor;
    }
    if (remainder >= divisor)
        ++quotient;
    return quotient;
}

MPUnit Physica::div2WordByFullWordR(MPUnit high, MPUnit low, MPUnit divisor) noexcept {
    assert(high < divisor);
    assert(divisor & MPUnitHighestBitMask);
    MPUnit quotient, quotient2;
    mulWordByWord(quotient, quotient2, high, getInverse(divisor));

    auto temp = quotient2;
    quotient2 += low;
    quotient += high + (temp > quotient2) + 1;
    auto remainder = low - quotient * divisor;
    if (remainder > quotient2)
        remainder += divisor;
    if (remainder >= divisor)
        remainder -= divisor;
    return remainder;
}

MPUnit Physica::divArrByFullArrWith1Word(const MPUnit* __restrict dividend, const MPUnit* __restrict divisor, size_t len) noexcept {
    assert(len >= 1);
    MPUnit q = dividend[len] >= divisor[len - 1] ? MPUnitMax : div2WordByFullWordQ(dividend[len], dividend[len - 1], divisor[len - 1]);
    if (len == 1) // May be ask len > 1 to avoid branches.
        return q;
    MPUnit temp[2]{dividend[len - 1], dividend[len]};
    /* Calculate temp - q * divisor[len - 1] */ {
        MPUnit temp_1[2];
        mulWordByWord(temp_1[1], temp_1[0], q, divisor[len - 1]);
        subArrByArrEq(temp, temp_1, 2);
    }
    MPUnit q_divisor_high, q_divisor_low;
    mulWordByWord(q_divisor_high, q_divisor_low, q, divisor[len - 2]);
    sub2WordByWord(q_divisor_high, q_divisor_low, dividend[len - 2]);
    if (temp[1] == 0 && (temp[0] < q_divisor_high || (temp[0] == q_divisor_high && q_divisor_low == 0)))
        --q;

    auto* n = HostAllocator<MPUnit>{}.allocate(len + 1);
    n[len] = mulArrByWord(n, divisor, len, q);
    // Judge dividend < n or not. If dividend < n, we have to do carry.
    bool carry = true;
    for (size_t i = len; i <= len; --i) { // i <= len makes use of overflow
        if (dividend[i] > n[i]) {
            carry = false;
            break;
        }
    }
    delete[] n;
    return carry ? (q - 1) : q;
}

std::unique_ptr<MPUnit[]> Physica::byteLeftShift(const MPUnit* __restrict byte, unsigned int length, unsigned int shift) noexcept {
    const unsigned int quotient = shift / MPUnitWidth;
    auto result = std::unique_ptr<MPUnit[]>(HostAllocator<MPUnit>{}.allocate(length));
    if (quotient < length) {
        if (quotient != 0) {
            memcpy(result.get() + quotient, byte, (length - quotient) * sizeof(MPUnit));
            memset(result.get(), 0, quotient * sizeof(MPUnit));
            shift -= quotient * MPUnitWidth;
        }
        if (shift != 0) {
            MPUnit carry = 0, temp;
            for (unsigned int i = quotient; i < length - 1; ++i) {
                temp = result[i] >> (MPUnitWidth - shift);
                result[i] <<= shift;
                result[i] |= carry;
                carry = temp;
            }
            result[length - 1] <<= shift;
            result[length - 1] |= carry;
        }
    }
    return result;
}

std::unique_ptr<MPUnit[]> Physica::byteRightShift(const MPUnit* __restrict byte, size_t length, size_t shift) noexcept {
    const size_t quotient = shift / MPUnitWidth;
    auto result = std::unique_ptr<MPUnit[]>(HostAllocator<MPUnit>{}.allocate(length));
    if (quotient < length) {
        auto bufferSize = length - quotient;
        if (quotient != 0) {
            memcpy(result.get(), byte + quotient, (length - quotient) * sizeof(MPUnit));
            memset(result.get() + bufferSize, 0, quotient * sizeof(MPUnit));
            shift -= quotient * MPUnitWidth;
        }
        if (shift != 0) {
            MPUnit carry = 0, temp;
            for (size_t i = bufferSize - 1; i > 0; --i) {
                temp = result[i] << (MPUnitWidth - shift);
                result[i] >>= shift;
                result[i] |= carry;
                carry = temp;
            }
            result[0] >>= shift;
            result[0] |= carry;
        }
    }
    return result;
}

void Physica::byteLeftShiftEq(MPUnit* __restrict byte, unsigned int length, unsigned int shift) noexcept {
    const unsigned int quotient = shift / MPUnitWidth;
    if (quotient < length) {
        if (quotient != 0) {
            memmove(byte + quotient, byte, (length - quotient) * sizeof(MPUnit));
            memset(byte, 0, quotient * sizeof(MPUnit));
            shift -= quotient * MPUnitWidth;
        }
        if (shift != 0) {
            MPUnit carry = 0, temp;
            for (unsigned int i = quotient; i < length - 1; ++i) {
                temp = byte[i] >> (MPUnitWidth - shift);
                byte[i] <<= shift;
                byte[i] |= carry;
                carry = temp;
            }
            byte[length - 1] <<= shift;
            byte[length - 1] |= carry;
        }
    }
}

void Physica::byteRightShiftEq(MPUnit* __restrict byte, size_t length, size_t shift) noexcept {
    const size_t quotient = shift / MPUnitWidth;
    if (quotient < length) {
        size_t bufferSize = length - quotient;
        if (quotient != 0) {
            memmove(byte, byte + quotient, bufferSize * sizeof(MPUnit));
            memset(byte + bufferSize, 0, quotient * sizeof(MPUnit));
            shift -= quotient * MPUnitWidth;
        }
        if (shift != 0) {
            MPUnit carry = 0, temp;
            for (size_t i = bufferSize - 1; i > 0; --i) {
                temp = byte[i] << (MPUnitWidth - shift);
                byte[i] >>= shift;
                byte[i] |= carry;
                carry = temp;
            }
            byte[0] >>= shift;
            byte[0] |= carry;
        }
    }
}
