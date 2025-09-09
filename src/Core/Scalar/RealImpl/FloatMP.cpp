/*
 * Copyright 2019-2025 Weibo He.
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
#include <iostream>
#include "Physica/Core/Scalar/Real.h"
#include "Physica/Core/Scalar/RealImpl/FloatMPImpl/AddBasic.h"
#include "Physica/Core/Scalar/RealImpl/FloatMPImpl/DivBasic.h"
#include "Physica/Core/Scalar/RealImpl/FloatMPImpl/ArraySupport.h"

using namespace Physica;

Real<FloatMP>::Real() : Real(GlobalPrecision, 0) {}

Real<FloatMP>::Real(int length_, int power_)
        : byte(HostAllocator<MPUnit>{}.allocate(std::abs(length_)))
        , length(length_)
        , power(power_) {
    assert(length != INT_MIN && "Length of scalar must not equal to INT_MIN or -length will make no sense");
}

Real<FloatMP>::Real(std::initializer_list<MPUnit> bytes_, int length_, int power_) : Real(length_, power_) {
    assert(bytes_.size() == std::abs(length_));
    auto *p = byte;
    for (auto elem : bytes_){
        *p = elem;
        p += 1;
    }
}

Real<FloatMP>::Real(SignedMPUnit x)
        : byte(HostAllocator<MPUnit>{}.allocate(1))
        , length(x > 0 ? 1 : -1)
        , power(0) {
    byte[0] = x > 0 ? x : -x;
}

Real<FloatMP>::Real(double d) noexcept {
    if (d == 0) {
        byte = HostAllocator<MPUnit>{}.allocate(1);
        length = 1;
        byte[0] = power = 0;
        return;
    }
    const double_extract extract{d};
    const int quotient = static_cast<int>(extract.exp) - 1023;
    power = quotient / int(MPUnitWidth);
    power -= quotient < 0; // Have power * MPUnitWidth < extract.exp, so that remainder > 0.
    unsigned int remainder = quotient - power * MPUnitWidth;
    if constexpr (PhysicaWordSize == 64) {
        if (remainder < 52) {
            length = 2;
            byte = HostAllocator<MPUnit>{}.allocate(length);
            // Hidden bit
            byte[1] = 1;
            byte[1] <<= remainder;
            if (remainder <= 20) {
                byte[1] += static_cast<MPUnit>(extract.high) >> (20 - remainder);
                byte[0] = static_cast<MPUnit>(extract.high) << (44 + remainder);
                byte[0] += static_cast<MPUnit>(extract.low) << (32 - (20 - remainder));
            }
            else {
                byte[1] += static_cast<MPUnit>(extract.high) << (remainder - 20);
                byte[1] += static_cast<MPUnit>(extract.low) >> (32 - (remainder - 20));
                byte[0] = static_cast<MPUnit>(extract.low) << (32 + (remainder - 20));
            }
        }
        else {
            length = 1;
            byte = HostAllocator<MPUnit>{}.allocate(1);
            // Hidden bit
            byte[0] = 1;
            byte[0] <<= 20U;
            byte[0] += static_cast<MPUnit>(extract.high);
            byte[0] <<= 32U;
            byte[0] += static_cast<MPUnit>(extract.low);
            byte[0] <<= remainder - 52;
        }
    }
    else {
        if (remainder < 20) {
            length = 3;
            byte = HostAllocator<MPUnit>{}.allocate(length);
            // Hidden bit
            byte[2] = 1;
            byte[2] <<= remainder;
            byte[2] += static_cast<MPUnit>(extract.high) >> (20 - remainder);
            byte[1] = static_cast<MPUnit>(extract.high) << (32 - (20 - remainder));
            byte[1] += static_cast<MPUnit>(extract.low) >> (20 - remainder);
            byte[0] = static_cast<MPUnit>(extract.low) << remainder;
        }
        else {
            length = 2;
            byte = HostAllocator<MPUnit>{}.allocate(length);
            // Hidden bit
            byte[1] = 1;
            byte[1] <<= remainder;
            byte[1] += static_cast<MPUnit>(extract.high) << (remainder - 20);
            byte[1] += static_cast<MPUnit>(extract.low) >> (32 - (remainder - 20));
            byte[0] = static_cast<MPUnit>(extract.low) << (remainder - 20);
        }
    }

    if (extract.sign)
        length = -length;
}

Real<FloatMP>::Real(const Integer& i)
        : byte(HostAllocator<MPUnit>{}.allocate(i.getSize()))
        , length(i.getLength())
        , power(i.getSize() - 1) {
    memcpy(byte, i.getByte(), getSize() * sizeof(MPUnit));
}

Real<FloatMP>::Real(const Rational& r) {
    Real<FloatMP> result = Real<FloatMP>(r.getNumerator()) / Real<FloatMP>(r.getDenominator());
    swap(result);
}
/**
 * Not accurate.
 */
Real<FloatMP>::Real(const char* s) : Real(strtod(s, nullptr)) {}
/**
 * Not accurate.
 */
Real<FloatMP>::Real(const wchar_t* s) {
    size_t size = wcslen(s);
    char* str = new char[size + 1];
    str[size] = '\0';
    for (size_t i = 0; i < size; ++i)
        str[i] = (char)s[i];
    Real<FloatMP> temp(str);
    byte = temp.byte;
    temp.byte = nullptr;
    length = temp.length;
    power = temp.power;
    delete[] str;
}

Real<FloatMP>::Real(const Real<FloatMP>& s)
        : byte(HostAllocator<MPUnit>{}.allocate(s.getSize()))
        , length(s.length)
        , power(s.power) {
    memcpy(byte, s.byte, getSize() * sizeof(MPUnit));
}

Real<FloatMP>::Real(Real<FloatMP>&& s) noexcept
        : byte(s.byte), length(s.length), power(s.power) {
    s.byte = nullptr;
}

Real<FloatMP>::~Real() {
    delete[] byte;
}

MPUnit Real<FloatMP>::operator[](unsigned int index) const {
    assert(index < static_cast<unsigned int>(getSize()));
    return byte[index];
}

Real<FloatMP>::operator double() const {
    if (isZero())
        return 0.0;
    return toDouble(byte, length, power);
}

Real<FloatMP> Real<FloatMP>::operator+(const Real<FloatMP>& s) const {
    auto result = add(*this, s);
    cutLength(result);
    return result;
}

Real<FloatMP> Real<FloatMP>::operator-(const Real<FloatMP>& s) const {
    auto result = sub(*this, s);
    cutLength(result);
    return result;
}

Real<FloatMP> Real<FloatMP>::operator*(const Real<FloatMP>& s) const {
    auto result = mul(*this, s);
    cutLength(result);
    return result;
}

Real<FloatMP> Real<FloatMP>::operator/(const Real<FloatMP>& s) const {
    return div(*this, s);
}

Real<FloatMP> Real<FloatMP>::operator<<(int bits) const {
    if (bits == 0)
        return Real(*this);
    if (bits < 0)
        return *this >> -bits;
    const int size = getSize();
    const int quotient = bits / MPUnitWidth; // NOLINT: quotient < INT_MAX
    const int remainder = bits - quotient * MPUnitWidth;
    // If remainder = 0, we must return directly because shifting a MPUnit for MPUnitWidth bits is a undefined behavior.
    if (remainder == 0) {
        Real copy(*this);
        copy.power += quotient;
        return copy;
    }

    const bool carry = std::countl_zero(byte[size - 1]) < remainder;
    Real result(length >= 0 ? (size + carry) : -(size + carry), power + quotient + carry);
    result.byte[0] = 0;
    const int size_1 = size - 1;
    for (int i = 0; i < size_1; ++i) {
        result.byte[i] |= byte[i] << remainder;
        result.byte[i + 1] = byte[i] >> (MPUnitWidth - remainder);
    }
    result.byte[size_1] |= byte[size_1] << remainder;
    if (carry)
        result.byte[size] = byte[size_1] >> (MPUnitWidth - remainder);
    return result;
}

Real<FloatMP> Real<FloatMP>::operator>>(int bits) const {
    if (bits == 0)
        return Real(*this);
    if (bits < 0)
        return *this << -bits;
    const int size = getSize();
    const int quotient = bits / MPUnitWidth; // NOLINT: quotient < INT_MAX
    const unsigned int remainder = bits - quotient * MPUnitWidth;
    // If remainder = 0, we must return directly because shifting a MPUnit for MPUnitWidth bits is a undefined behavior.
    if (remainder == 0) {
        Real copy(*this);
        copy.power -= quotient;
        return copy;
    }

    const bool carry = (std::countl_zero(byte[size - 1]) + remainder) < MPUnitWidth;
    Real result(length >= 0 ? (size + carry) : -(size + carry), power - quotient + carry - 1);
    if (carry)
        result.byte[size] = byte[size - 1] >> remainder;

    for (int i = size - 1; i > 0; --i) {
        result.byte[i] = byte[i] << (MPUnitWidth - remainder);
        result.byte[i] |= byte[i - 1] >> remainder;
    }
    result.byte[0] = byte[0] << (MPUnitWidth - remainder);
    return result;
}

Real<FloatMP> Real<FloatMP>::operator-() const {
    Real result(-length, power);
    memcpy(result.byte, byte, getSize() * sizeof(MPUnit));
    return result;
}
/*!
 * Optimize:
 * Is subtract faster than comparing the elements?
 */
bool Real<FloatMP>::operator>(const This& other) const {
    // Judge from sign.
    const bool positive = isPositive();
    if (positive) {
        if (!other.isPositive())
            return true;
    }
    else if (isZero())
        return other.isNegative();
    else {
        if (!other.isNegative())
            return false;
    }
    // If we cannot get a result, judge from power
    bool result;
    if (getPower() > other.getPower())
        result = positive;
    else if (getPower() < other.getPower())
        result = !positive;
    else {
        // The only method left.
        // Optimize: We have confirmed that both scalars have the same sign and power, possible make use them to get better performance.
        Real<FloatMP> subtract = *this - other;
        result = subtract.isPositive();
    }
    return result;
}
/*!
 * Optimize:
 * Is subtract faster than comparing the elements?
 */
bool Real<FloatMP>::operator<(const This& other) const {
    // Judge from sign.
    const bool positive = isPositive();
    if (positive) {
        if (!other.isPositive())
            return false;
    }
    else if (isZero())
        return other.isPositive();
    else {
        if (!other.isNegative())
            return true;
    }
    // If we cannot get a result, judge from power
    bool result;
    if (getPower() > other.getPower())
        result = !positive;
    else if (getPower() < other.getPower())
        result = positive;
    else {
        // The only method left.
        // Optimize: We have confirmed that both scalars have the same sign and power, possible make use them to get better performance.
        Real<FloatMP> subtract = *this - other;
        result = subtract.isNegative();
    }
    return result;
}

bool Real<FloatMP>::operator==(const This& other) const {
    return getPower() == other.getPower()
        // Here length may not equal n.length because we define numbers like 1.0 and 1.00 are equal.
        && ((getLength() ^ other.getLength()) >= 0) // NOLINT
        // Optimize: We have confirmed that both scalars have the same sign and power, possible make use them to get better performance.
        && Real<FloatMP>(*this - other).isZero();
}

void Real<FloatMP>::swap(This& __restrict s) noexcept {
    assert(this != &s && "[Error]: Self swap is likely a bug");
    std::swap(byte, s.byte);
    std::swap(length, s.length);
    std::swap(power, s.power);
}

void Real<FloatMP>::dump() const noexcept {
    for (int i = 0; i < std::abs(length); ++i)
        std::cout << byte[i] << ' ';
    std::cout << std::format("\n{} {}\n", length, power);
}
/**
 * Returns true if s1 and s2 has the same sign. Both s1 and s2 do not equal to zero.
 * This function provide a quick sign check compare to using isPositive() and isNegative().
 */
bool Real<FloatMP>::matchSign(const Real<FloatMP>& s1, const Real<FloatMP>& s2) noexcept {
    assert(!s1.isZero() && !s2.isZero());
    return (s1.length ^ s2.length) >= 0; // NOLINT Bitwise operator between two signed integer is intended.
}

double Real<FloatMP>::toDouble(MPUnit* __restrict byte, int length, int power) noexcept {
    double_extract extract{0};
    extract.sign = length < 0;

    const auto size = std::abs(length);
    const auto zeroCount = std::countl_zero(byte[size - 1]); // Optimize: (size - 1) and (size > 1) is used several times
    // Using long to avoid overflow.
    const long exp = power * int(MPUnitWidth) + static_cast<long>(MPUnitWidth - zeroCount) - 1 + 1023;
    if (exp >= 2047) {
        extract.high = extract.low = 0;
        extract.exp = 2047;
        return extract.value;
    }

    if (exp <= 0)
        return 0.0;
    extract.exp = exp;

    auto temp = byte[size - 1] << (zeroCount + 1);
    if constexpr (PhysicaWordSize == 64) {
        extract.high = temp >> 44U;
        if (zeroCount <= 11) {
            extract.low = temp << 20U >> 32U;
        }
        else {
            if (zeroCount <= 44 - 1) {
                extract.low = temp << 20U >> 32U;
                if (size > 1)
                    extract.low += byte[size - 2] >> (64 - (32 - (64 - 20 - zeroCount - 1)));
            }
            else {
                if (size > 1) {
                    extract.high += byte[size - 2] >> (64 - (20 - (64 - zeroCount - 1)));
                    extract.low = byte[size - 2] << (20 - (64 - zeroCount - 1)) >> 32U;
                }
            }
        }
    }
    else {
        extract.high = temp >> 12U;
        if (zeroCount <= 11) {
            extract.low = temp << 20U;
            if (size > 1)
                extract.low = byte[size - 1] >> (32 - 20 - zeroCount - 1);
        }
        else {
            if (size > 1) {
                extract.high += byte[size - 1] >> (32 - (zeroCount + 1 - 12));
                extract.low = byte[size - 1] << (zeroCount + 1 - 12);
            }
            if (size > 2)
                extract.low += byte[size - 2] >> (32 - (zeroCount + 1 - 12));
        }
    }
    return extract.value;
}
/**
 * Cut zeros from the beginning.
 */
void Real<FloatMP>::cutZero() {
    const int size = getSize();
    int id = size - 1;
    while (byte[id] == 0 && id > 0)
        --id;
    ++id;

    if (id != size) {
        int shorten = size - id;
        byte = reinterpret_cast<MPUnit*>(realloc(byte, id * sizeof(MPUnit)));
        length = length > 0 ? id : -id;
        auto temp = power;
        power = byte[id - 1] != 0 ? (temp - shorten) : 0;
    }
}

Real<FloatMP> Real<FloatMP>::add(const Real<FloatMP>& s1, const Real<FloatMP>& s2) noexcept {
    if (s1.isZero()) [[unlikely]]
        return s2;
    if (s2.isZero()) [[unlikely]]
        return s1;

    if (!matchSign(s1, s2)) {
        if (s1.length > 0) {
            Real shallow_copy(const_cast<MPUnit*>(s2.byte), -s2.length, s2.power);
            auto result = sub(s1, shallow_copy);
            shallow_copy.byte = nullptr;
            return result;
        }
        else {
            Real shallow_copy(const_cast<MPUnit*>(s1.byte), -s1.length, s1.power);
            auto result = sub(s2, shallow_copy);
            shallow_copy.byte = nullptr;
            return result;
        }
    }

    const bool flag = s1.power > s2.power;
    const Real& big = flag ? s1 : s2;
    const Real& small = flag ? s2 : s1;
    const int bigSize = big.getSize();
    const int smallSize = small.getSize();
    // Estimate the ed of result first, will calculate it accurately later.
    int length = (big.power + std::max(bigSize - big.power, smallSize - small.power));
    length = length > GlobalPrecision
                    ? GlobalPrecision
                    : length;
    auto* byte = HostAllocator<MPUnit>{}.allocate(length);
    /* Init byte */ {
        const auto copySize = bigSize > length ? length : bigSize;
        const auto clearSize = length - copySize;
        memset(byte, 0, clearSize * sizeof(MPUnit));
        memcpy(byte + clearSize, big.byte + bigSize - copySize, copySize * sizeof(MPUnit));
    }
    bool carry;
    /* Add and carry */ {
        const int dist = small.power - (big.power - GlobalPrecision);
        const int overlap = dist > smallSize ? smallSize
                                             : ((dist < 0) ? 0 : dist);
        carry = addArrWithArrEq(small.byte + smallSize - overlap, byte + length + small.power - big.power - overlap, overlap);

        int index = dist;
        while (carry != 0 && index < length) {
            MPUnit temp = byte[index] + 1;
            byte[index] = temp;
            carry = temp == 0;
            ++index;
        }
    }
    ///////////////////////////////////////Get byte, length and power//////////////////////////
    int power = big.power;
    if (carry) {
        ++length;
        ++power;
        byte = reinterpret_cast<MPUnit*>(realloc(byte, length * sizeof(MPUnit)));
        byte[length - 1] = 1;
    }
    return Real<FloatMP>(byte, big.length < 0 ? -length : length, power);
}

Real<FloatMP> Real<FloatMP>::sub(const Real<FloatMP>& s1, const Real<FloatMP>& s2) noexcept {
    if (s1.isZero()) [[unlikely]]
        return -s2;
    if (s2.isZero()) [[unlikely]]
        return s1;

    if (s1.length > 0) {
        if (s2.length < 0) {
            Real shallow_copy(const_cast<MPUnit*>(s2.byte), -s2.length, s2.power);
            Real result = add(s1, shallow_copy);
            shallow_copy.byte = nullptr;
            return result;
        }
        else {
            const Real* big;
            const Real* small;
            bool changeSign = false;
            if (s1.power > s2.power) {
                big = &s1;
                small = &s2;
            }
            else {
                changeSign = true;
                big = &s2;
                small = &s1;
            }
        redo:
            const int bigSize = big->getSize();
            const int smallSize = small->getSize();
            // Estimate the ed of result first, will calculate it accurately later.
            int length = (big->power + std::max(bigSize - big->power, smallSize - small->power));
            length = length > GlobalPrecision
                           ? GlobalPrecision
                           : length;
            auto* byte = HostAllocator<MPUnit>{}.allocate(length);
            /* Init byte */ {
                const auto copySize = bigSize > length ? length : bigSize;
                const auto clearSize = length - copySize;
                memset(byte, 0, clearSize * sizeof(MPUnit));
                memcpy(byte + clearSize, big->byte + bigSize - copySize, copySize * sizeof(MPUnit));
            }
            bool carry;
            MPUnit a;
            /* Sub and carry */ {
                // usableSmall is the part whose sub result will fall in GlobalPrecision.
                int usableSmall = small->power - (big->power - GlobalPrecision);
                a = usableSmall < 0;
                usableSmall = usableSmall > smallSize
                                    ? smallSize
                                    : (a ? 0 : usableSmall);
                int smallEnd = length + small->power - big->power;
                carry = subArrByArrEq(byte + smallEnd - usableSmall, small->byte + smallSize - usableSmall, usableSmall);
                // usableSmall is also the index which we should carry to.
                MPUnit temp1, temp2;
                while (carry != 0 && smallEnd < length) {
                    temp1 = byte[smallEnd];
                    temp2 = temp1 - 1;
                    byte[smallEnd] = temp2;
                    carry = temp1 < temp2;
                    ++smallEnd;
                }
            }

            if (carry) {
                auto temp = big;
                big = small;
                small = temp;
                changeSign = !changeSign;
                delete[] byte;
                goto redo;
            }
            Real<FloatMP> result(byte, changeSign ? -length : length, big->power);
            result.cutZero();
            return result;
        }
    }
    else {
        Real shallow_copy(const_cast<MPUnit*>(s1.byte), -s1.length, s1.power);
        if (s2.length > 0) {
            Real result = add(shallow_copy, s2);
            result.toOpposite();
            shallow_copy.byte = nullptr;
            return result;
        }
        else {
            Real shallow_copy_1(const_cast<MPUnit*>(s2.byte), -s2.length, s2.power);
            Real result = sub(shallow_copy_1, shallow_copy);
            shallow_copy.byte = shallow_copy_1.byte = nullptr;
            return result;
        }
    }
}
// Optimize: length may be too long and it is unnecessary, cut it and consider the accuracy.
Real<FloatMP> Real<FloatMP>::mul(const Real<FloatMP>& s1, const Real<FloatMP>& s2) noexcept {
    if (s1.isZero() || s2.isZero()) [[unlikely]]
        return Real<FloatMP>(0);
    if (s1 == BasicConst::getInstance()._1) [[unlikely]]
        return s2;
    if (s2 == BasicConst::getInstance()._1) [[unlikely]]
        return s1;
    const int size1 = s1.getSize();
    const int size2 = s2.getSize();
    // Estimate the ed of result first. we will calculate it accurately later.
    auto length = size1 + size2;
    auto* byte = HostAllocator<MPUnit>{}.allocate(length);
    for (int i = 0; i < size2; ++i)
        byte[i] = 0;
    for (int i = 0; i < size1; ++i)
        byte[i + size2] = mulAddArrByWord(byte + i, s2.byte, size2, s1.byte[i]);
    ///////////////////////////////////////Get byte, length and power//////////////////////////;
    int power = s1.power + s2.power + 1;
    if (byte[length - 1] == 0) {
        --length;
        --power;
        byte = reinterpret_cast<MPUnit*>(realloc(byte, length * sizeof(MPUnit)));
    }
    ////////////////////////////////////Out put////////////////////////////////////////
    return Real<FloatMP>(byte, matchSign(s1, s2) ? length : -length, power);
}

Real<FloatMP> Real<FloatMP>::div(const Real<FloatMP>& s1, const Real<FloatMP>& s2) noexcept {
    assert(!s2.isZero() && "[Error]: Divide by zero");
    if (!s1.isZero()) [[likely]] {
        if (s2 != BasicConst::getInstance()._1) {
            const auto s1_size = s1.getSize(), s2_size = s2.getSize();
            // Add one to arr1_length to avoid precision loss during right shift.
            auto arr1_len = std::max(s1_size, s2_size) + 1;
            auto s1_blank = arr1_len - s1_size;
            auto arr1 = HostAllocator<MPUnit>{}.allocate(arr1_len);
            memcpy(arr1 + s1_blank, s1.byte, s1_size * sizeof(MPUnit));
            memset(arr1, 0, s1_blank * sizeof(MPUnit));
            // Size of arr2 is arranged 1 less than arr1.
            auto arr2_len = arr1_len - 1;
            auto s2_blank = arr2_len - s2_size;
            auto arr2 = HostAllocator<MPUnit>{}.allocate(arr2_len);
            memcpy(arr2 + s2_blank, s2.byte, s2_size * sizeof(MPUnit));
            memset(arr2, 0, s2_blank * sizeof(MPUnit));
            /*
             * We shift s1 and s2, making the less highest bit of s1 is set and the highest bit of s2 is set
             * to meet the acquirement of the function divArrByFullArrWith1Word().
             */
            const int s1_shift = static_cast<int>(std::countl_zero(s1.byte[s1_size - 1])) - 1;
            if (s1_shift > 0)
                byteLeftShiftEq(arr1, arr1_len, s1_shift);
            else
                byteRightShiftEq(arr1, arr1_len, -s1_shift);
            const int s2_shift = static_cast<int>(std::countl_zero(s2.byte[s2_size - 1]));
            byteLeftShiftEq(arr2, arr2_len, s2_shift);
            ////////////////////////////////Calculate cursory first//////////////////////////////////////
            // Estimate the length of result.
            int length = GlobalPrecision;
            auto* byte = HostAllocator<MPUnit>{}.allocate(length);

            for (int i = length - 1; i >= 0; --i) {
                byte[i] = divArrByFullArrWith1Word(arr1, arr2, arr2_len);
                arr1[arr2_len] -= mulSubArrByWord(arr1, arr2, arr2_len, byte[i]);
                byteLeftShiftEq(arr1, arr1_len, MPUnitWidth);
            }
            delete[] arr1;
            delete[] arr2;
            ////////////////////////////////////Out put////////////////////////////////////////
            return Real<FloatMP>(byte, matchSign(s1, s2) ? length : -length, s1.getPower() - s2.getPower() - 1) >> (s1_shift - s2_shift);
        }
        else
            return static_cast<const Real<FloatMP>&>(s1);
    }
    else
        return Real<FloatMP>(static_cast<SignedMPUnit>(0));
}
/**
 * If the length of new array is larger than GlobalPrecision, it will be set to GlobalPrecision.
 * Return true if array is cut.
 */
bool Real<FloatMP>::cutLength(Real<FloatMP>& s) {
    bool result = false;
    int size = s.getSize();

    if (size > GlobalPrecision) {
        result = true;
        int cutFrom = size - GlobalPrecision;
        auto new_byte = HostAllocator<MPUnit>{}.allocate(GlobalPrecision);
        memcpy(new_byte, s.byte + cutFrom, GlobalPrecision * sizeof(MPUnit));
        delete[] s.byte;
        s.byte = new_byte;
        s.length = s.length > 0 ? GlobalPrecision : -GlobalPrecision;
    }
    return result;
}

namespace Physica {
    /*!
     * return true if the abstract value of s1 is larger or equal than the abstract value of s2.
     * return false if the abstract value of s1 is smaller to the abstract value of s2.
     *
     * Optimize:
     * Is subtract faster than comparing the elements?
     */
    bool absCompare(const Real<FloatMP>& s1, const Real<FloatMP>& s2) {
        if (s1.isZero() || s2.isZero())
            return true;
        if (s1.getPower() > s2.getPower())
            return true;
        if (s1.getPower() < s2.getPower())
            return false;
        const Real<FloatMP>*longer, *shorter;
        int longer_length, shorter_length;
        /* Compare length */ {
            const auto n1_length = s1.getLength(), n2_length = s2.getLength();
            bool b = n1_length > n2_length;
            longer = b ? &s1 : &s2;
            shorter = b ? &s2 : &s1;
            longer_length = b ? n1_length : n2_length;
            shorter_length = b ? n2_length : n1_length;
        }
        --longer_length;
        --shorter_length;
        for (; shorter_length >= 0; --shorter_length, --longer_length) {
            if ((*longer)[longer_length] > (*shorter)[shorter_length])
                return longer == &s1;
            if ((*longer)[longer_length] < (*shorter)[shorter_length])
                return shorter == &s1;
        }
        return longer_length == shorter_length || longer == &s1;
    }

    Real<FloatMP>& operator++(Real<FloatMP>& s) {
        s += BasicConst::getInstance()._1;
        return s;
    }

    Real<FloatMP>& operator--(Real<FloatMP>& s) {
        s -= BasicConst::getInstance()._1;
        return s;
    }

    Real<FloatMP> operator++(Real<FloatMP>& s, int) {
        Real<FloatMP> temp(s);
        s += BasicConst::getInstance()._1;
        return temp;
    }

    Real<FloatMP> operator--(Real<FloatMP>& s, int) {
        Real<FloatMP> temp(s);
        s -= BasicConst::getInstance()._1;
        return temp;
    }
}
