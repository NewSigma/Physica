/*
 * Copyright 2019 WeiBo He.
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
#include "Physica/Core/MultiPrecition/Scalar.h"

namespace Physica::Core {
    //////////////////////////////////MultiPrecision-WithoutError///////////////////////////////////
    Scalar<MultiPrecision, false>::Scalar() noexcept
            : byte(nullptr), length(0), power(0) {}

    Scalar<MultiPrecision, false>::Scalar(int length, int power) noexcept
            : byte(reinterpret_cast<ScalarUnit*>(malloc(abs(length) * sizeof(ScalarUnit))))
            , length(length), power(power) {}

    Scalar<MultiPrecision, false>::Scalar(const Scalar<MultiPrecision, false>& s)
            : byte(reinterpret_cast<ScalarUnit*>(malloc(s.getSize() * sizeof(ScalarUnit))))
            , length(s.length), power(s.power) {
        memcpy(byte, s.byte, getSize() * sizeof(ScalarUnit));
    }

    Scalar<MultiPrecision, false>::Scalar(Scalar<MultiPrecision, false>&& s) noexcept
            : byte(s.byte), length(s.length), power(s.power) {
        s.byte = nullptr;
    }

    Scalar<MultiPrecision, false>::Scalar(int i) : Scalar(static_cast<SignedScalarUnit>(i)) {}

    Scalar<MultiPrecision, false>::Scalar(SignedScalarUnit unit)
            : byte(reinterpret_cast<ScalarUnit*>(malloc(sizeof(ScalarUnit))))
            , length(unit > 0 ? 1 : -1), power(0) {
        byte[0] = unit > 0 ? unit : -unit;
    }

    Scalar<MultiPrecision, false>::Scalar(double d) {
        if(d == 0) {
            byte = reinterpret_cast<ScalarUnit*>(malloc(sizeof(ScalarUnit)));
            length = 1;
            byte[0] = power = 0;
            return;
        }
        double_extract extract{d};
        auto quotient = static_cast<int>(extract.exp) - 1023;
        power = quotient / __WORDSIZE;
        //Have power * __WORDSIZE < extract.exp, so that remainder > 0.
        if(quotient < 0)
            --power;
        unsigned int remainder = quotient - power * __WORDSIZE;
    #ifdef PHYSICA_64BIT
        if(remainder < 52) {
            length = 2;
            byte = reinterpret_cast<ScalarUnit*>(malloc(length * sizeof(ScalarUnit)));
            //Hidden bit
            byte[1] = 1;
            byte[1] <<= remainder;
            if(remainder <= 20) {
                byte[1] += static_cast<ScalarUnit>(extract.high) >> (20 - remainder);
                byte[0] = static_cast<ScalarUnit>(extract.high) << (44 + remainder);
                byte[0] += static_cast<ScalarUnit>(extract.low) << (32 - (20 - remainder));
            }
            else {
                byte[1] += static_cast<ScalarUnit>(extract.high) << (remainder - 20);
                byte[1] += static_cast<ScalarUnit>(extract.low) >> (32 - (remainder - 20));
                byte[0] = static_cast<ScalarUnit>(extract.low) << (32 + (remainder - 20));
            }
        }
        else {
            length = 1;
            byte = reinterpret_cast<ScalarUnit*>(malloc(sizeof(ScalarUnit)));
            //Hidden bit
            byte[0] = 1;
            byte[0] <<= 20U;
            byte[0] += static_cast<ScalarUnit>(extract.high);
            byte[0] <<= 32U;
            byte[0] += static_cast<ScalarUnit>(extract.low);
            byte[0] <<= remainder - 52;
        }
    #endif
    #ifdef PHYSICA_32BIT
        if(remainder < 20) {
            length = 3;
            byte = reinterpret_cast<ScalarUnit*>(malloc(length * sizeof(ScalarUnit)));
            //Hidden bit
            byte[2] = 1;
            byte[2] <<= remainder;
            byte[2] += static_cast<ScalarUnit>(extract.high) >> (20 - remainder);
            byte[1] = static_cast<ScalarUnit>(extract.high) << (32 - (20 - remainder));
            byte[1] +=  static_cast<ScalarUnit>(extract.low) >> (20 - remainder);
            byte[0] = static_cast<ScalarUnit>(extract.low) << remainder;
        }
        else {
            length = 2;
            byte = reinterpret_cast<ScalarUnit*>(malloc(length * sizeof(ScalarUnit)));
            //Hidden bit
            byte[1] = 1;
            byte[1] <<= remainder;
            byte[1] += static_cast<ScalarUnit>(extract.high) << (remainder - 20);
            byte[1] += static_cast<ScalarUnit>(extract.low) >> (32 - (remainder - 20));
            byte[0] = static_cast<ScalarUnit>(extract.low) << (remainder - 20);
        }
    #endif
        if(extract.sign)
            length = -length;
    }
    /*!
     * Not accurate.
     */
    Scalar<MultiPrecision, false>::Scalar(const char* s) : Scalar(strtod(s, nullptr)) {}
    /*!
     * Not accurate.
     */
    Scalar<MultiPrecision, false>::Scalar(const wchar_t* s) {
        size_t size = wcslen(s);
        char str[size + 1];
        str[size] = '\0';
        for(int i = 0; i < size; ++i)
            str[i] = (char)s[i];
        Scalar temp(str);
        byte = temp.byte;
        temp.byte = nullptr;
        length = temp.length;
        power = temp.power;
    }

    Scalar<MultiPrecision, false>::~Scalar() { free(byte); }

    Scalar<MultiPrecision, false>& Scalar<MultiPrecision, false>::operator= (
            const Scalar<MultiPrecision, false>& s) {
        if(this == &s)
            return *this;
        length = s.length;
        int size = getSize();
        this->~Scalar();
        byte = reinterpret_cast<ScalarUnit*>(malloc(size * sizeof(ScalarUnit)));
        memcpy(byte, s.byte, size * sizeof(ScalarUnit));
        power = s.power;
        return *this;
    }

    Scalar<MultiPrecision, false>& Scalar<MultiPrecision, false>::operator=(
            Scalar<MultiPrecision, false>&& s) noexcept {
        this->~Scalar();
        byte = s.byte;
        s.byte = nullptr;
        length = s.length;
        power = s.power;
        return *this;
    }

    Scalar<MultiPrecision, false>::operator double() const {
        if(isZero())
            return 0.0;
        double_extract extract{0};
        extract.sign = length < 0;

        const auto zeroCount = countLeadingZeros(byte[getSize() - 1]);
        //Using long to avoid overflow.
        const long exp = power * __WORDSIZE + static_cast<long>(ScalarUnitWidth - zeroCount) - 1 + 1023;
        if(exp >= 2047) {
            extract.high = extract.low = 0;
            extract.exp = 2047;
            return extract.value;
        }
        if(exp <= 0) {
            return 0.0;
        }
        extract.exp = exp;

        auto size = getSize();
        auto temp = byte[size - 1] << (zeroCount + 1);
    #ifdef PHYSICA_64BIT
        extract.high = temp >> 44U;
        if(zeroCount <= 11) {
            extract.low = temp << 20U >> 32U;
        }
        else {
            if(zeroCount <= 44 - 1) {
                extract.low = temp << 20U >> 32U;
                if(size > 1)
                    extract.low += byte[size - 2] >> (64 - (32 - (64 - 20 - zeroCount - 1)));
            }
            else {
                if(size > 1) {
                    extract.high += byte[size - 2] >> (64 - (20 - (64 - zeroCount - 1)));
                    extract.low = byte[size - 2] << (20 - (64 - zeroCount - 1)) >> 32U;
                }
            }
        }
    #endif
    #ifdef PHYSICA_32BIT
        extract.high = temp >> 12U;
        if(zeroCount <= 11) {
            extract.low = temp << 20U;
            if(size > 1)
                extract.low = byte[size - 1] >> (32 - 20 - zeroCount - 1);
        }
        else {
            if(size > 1) {
                extract.high += byte[size - 1] >> (32 - (zeroCount + 1 - 12));
                extract.low = byte[size - 1] << (zeroCount + 1 - 12);
            }
            if(size > 2)
                extract.low += byte[size - 2] >> (32 - (zeroCount + 1 - 12));
        }
    #endif
        return extract.value;
    }

    Scalar<MultiPrecision, false> Scalar<MultiPrecision, false>::operator+(
            const Scalar<MultiPrecision, false>& s) const {
        auto result = add<false>(*this, s);
        cutLength(result);
        return result;
    }

    Scalar<MultiPrecision, false> Scalar<MultiPrecision, false>::operator-(
            const Scalar<MultiPrecision, false>& s) const {
        auto result = sub<false>(*this, s);
        cutLength(result);
        return result;
    }

    Scalar<MultiPrecision, false> Scalar<MultiPrecision, false>::operator*(
            const Scalar<MultiPrecision, false>& s) const {
        auto result = mul<false>(*this, s);
        cutLength(result);
        return result;
    }

    Scalar<MultiPrecision, false> Scalar<MultiPrecision, false>::operator/(
            const Scalar<MultiPrecision, false>& s) const {
        return div<false>(*this, s);
    }

    Scalar<MultiPrecision, true> Scalar<MultiPrecision, false>::operator+(
            const Scalar<MultiPrecision, true>& s) const {
        auto result = add<true>(*this, s);
        cutLength(result);
        if(s.getA() != 0)
            result.applyError(s.getAccuracy());
        return result;
    }

    Scalar<MultiPrecision, true> Scalar<MultiPrecision, false>::operator-(
            const Scalar<MultiPrecision, true>& s) const {
        auto result = sub<true>(*this, s);
        cutLength(result);
        if(s.getA() != 0)
            result.applyError(s.getAccuracy());
        return result;
    }

    Scalar<MultiPrecision, true> Scalar<MultiPrecision, false>::operator*(
            const Scalar<MultiPrecision, true>& s) const {
        auto result = mul<true>(*this, s);
        cutLength(result);
        if(s.getA() != 0)
            result.applyError(*this * s.getAccuracy());
        return result;
    }

    Scalar<MultiPrecision, true> Scalar<MultiPrecision, false>::operator/(
            const Scalar<MultiPrecision, true>& s) const {
        auto result = div<true>(*this, s);
        if(s.getA() != 0) {
            const Scalar<MultiPrecision, false>& casted = s;
            Scalar s_a = s.getAccuracy();
            Scalar temp_1 = *this * s_a;
            Scalar temp_2 = casted * (casted - s_a);
            result.applyError(temp_1 / temp_2);
        }
        return result;
    }

    Scalar<MultiPrecision, false> Scalar<MultiPrecision, false>::operator<<(int bits) const {
        if(bits == 0)
            return Scalar(*this);
        if(bits < 0)
            return *this >> -bits;
        const int size = getSize();
        const int quotient = bits / ScalarUnitWidth; //NOLINT: quotient < INT_MAX
        const unsigned int remainder = bits - quotient * ScalarUnitWidth;
        //If remainder = 0, we must return directly because shifting a ScalarUnit for ScalarUnitWidth bits is a undefined behavior.
        if(remainder == 0) {
            Scalar copy(*this);
            copy.power += quotient;
            return copy;
        }

        const bool carry = countLeadingZeros(byte[size - 1]) < remainder;
        Scalar result(length >= 0 ? (size + carry) : -(size + carry), power + quotient + carry);
        result[0] = 0;
        const int size_1 = size - 1;
        for(int i = 0; i < size_1; ++i) {
            result[i] |= byte[i] << remainder;
            result[i + 1] = byte[i] >> (ScalarUnitWidth - remainder);
        }
        result[size_1] |= byte[size_1] << remainder;
        if(carry)
            result[size] = byte[size_1] >> (ScalarUnitWidth - remainder);
        return result;
    }

    Scalar<MultiPrecision, false> Scalar<MultiPrecision, false>::operator>>(int bits) const {
        if(bits == 0)
            return Scalar(*this);
        if(bits < 0)
            return *this << -bits;
        const int size = getSize();
        const int quotient = bits / ScalarUnitWidth; //NOLINT: quotient < INT_MAX
        const unsigned int remainder = bits - quotient * ScalarUnitWidth;
        //If remainder = 0, we must return directly because shifting a ScalarUnit for ScalarUnitWidth bits is a undefined behavior.
        if(remainder == 0) {
            Scalar copy(*this);
            copy.power -= quotient;
            return copy;
        }

        const bool carry = (countLeadingZeros(byte[size - 1]) + remainder) < ScalarUnitWidth;
        Scalar result(length >= 0 ? (size + carry) : -(size + carry), power - quotient + carry - 1);
        if(carry)
            result[size] = byte[size - 1] >> remainder;

        for(int i = size - 1; i > 0; --i) {
            result[i] = byte[i] << (ScalarUnitWidth - remainder);
            result[i] |= byte[i - 1] >> remainder;
        }
        result[0] = byte[0] << (ScalarUnitWidth - remainder);
        return result;
    }

    Scalar<MultiPrecision, false> Scalar<MultiPrecision, false>::operator-() const {
        Scalar result(-length, power);
        memcpy(result.byte, byte, getSize() * sizeof(ScalarUnit));
        return result;
    }

    /*!
     * return true if the abstract value of s1 is larger or equal than the abstract value of s2.
     * return false if the abstract value of s1 is smaller to the abstract value of s2.
     *
     * Optimize:
     * Is subtract faster than comparing the elements?
     */
    bool absCompare(const Scalar<MultiPrecision, false>& s1, const Scalar<MultiPrecision, false>& s2) {
        if(s1.isZero() || s2.isZero())
            return true;
        if(s1.getPower() > s2.getPower())
            return true;
        if(s1.getPower() < s2.getPower())
            return false;
        const Scalar<MultiPrecision, false>* longer, *shorter;
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
        for(; shorter_length >= 0; --shorter_length, --longer_length) {
            if((*longer)[longer_length] > (*shorter)[shorter_length])
                return longer == &s1;
            if((*longer)[longer_length] < (*shorter)[shorter_length])
                return shorter == &s1;
        }
        return longer_length == shorter_length || longer == &s1;
    }
    /*!
     * Optimize:
     * Is subtract faster than comparing the elements?
     */
    bool operator> (const Scalar<MultiPrecision, false>& s1, const Scalar<MultiPrecision, false>& s2) {
        //Judge from sign.
        const bool positive = s1.isPositive();
        if(positive) {
            if(!s2.isPositive())
                return true;
        }
        else if(s1.isZero())
            return s2.isNegative();
        else {
            if(!s2.isNegative())
                return false;
        }
        //If we cannot get a result, judge from power
        bool result;
        if(s1.getPower() > s2.getPower())
            result = positive;
        else if(s1.getPower() < s2.getPower())
            result = !positive;
        else {
            //The only method left.
            //Optimize: We have confirmed that s1, s2 have the same sign and power, possible make use them to get better performance.
            Scalar<MultiPrecision, false> subtract = s1 - s2;
            result = subtract.isPositive();
        }
        return result;
    }
    /*!
     * Optimize:
     * Is subtract faster than comparing the elements?
     */
    bool operator< (const Scalar<MultiPrecision, false>& s1, const Scalar<MultiPrecision, false>& s2) {
        //Judge from sign.
        const bool positive = s1.isPositive();
        if(positive) {
            if(!s2.isPositive())
                return false;
        }
        else if(s1.isZero())
            return s2.isPositive();
        else {
            if(!s2.isNegative())
                return true;
        }
        //If we cannot get a result, judge from power
        bool result;
        if(s1.getPower() > s2.getPower())
            result = !positive;
        else if(s1.getPower() < s2.getPower())
            result = positive;
        else {
            //The only method left.
            //Optimize: We have confirmed that s1, s2 have the same sign and power, possible make use them to get better performance.
            Scalar<MultiPrecision, false> subtract = s1 - s2;
            result = subtract.isNegative();
        }
        return result;
    }

    bool operator== (const Scalar<MultiPrecision, false>& s1, const Scalar<MultiPrecision, false>& s2) {
        return s1.getPower() == s2.getPower()
               //Here length may not equal n.length because we define numbers like 1.0 and 1.00 are equal.
               && ((s1.getLength() ^ s2.getLength()) >= 0) //NOLINT
               //Optimize: We have confirmed that s1, s2 have the same sign and power, possible make use them to get better performance.
               && (s1 - s2).isZero();
    }

    void Scalar<MultiPrecision, false>::swap(Scalar<MultiPrecision, false>& s) noexcept {
        std::swap(byte, s.byte);
        std::swap(length, s.length);
        std::swap(power, s.power);
    }
    /*!
     * Return true if s1 and s2 has the same sign. Both s1 and s2 do not equal to zero.
     * This function provide a quick sign check compare to using isPositive() and isNegative().
     */
    inline bool Scalar<MultiPrecision, false>::matchSign(
            const Scalar<MultiPrecision, false>& s1, const Scalar<MultiPrecision, false>& s2) {
        Q_ASSERT(!s1.isZero() && !s2.isZero());
        return (s1.length ^ s2.length) >= 0; //NOLINT Bitwise operator between two signed integer is intended.
    }
    ///////////////////////////////////MultiPrecision-WithError/////////////////////////////////////
    Scalar<MultiPrecision, true>::Scalar() noexcept : Scalar<MultiPrecision, false>(), a(0) {}

    Scalar<MultiPrecision, true>::Scalar(int length, int power, ScalarUnit a) noexcept
            : Scalar<MultiPrecision, false>(length, power), a(a) {}

    Scalar<MultiPrecision, true>::Scalar(const Scalar<MultiPrecision, true>& s) : Scalar<MultiPrecision, false>(s) {
        a = s.a;
    }

    Scalar<MultiPrecision, true>::Scalar(Scalar<MultiPrecision, true>&& s) noexcept
            : Scalar<MultiPrecision, false>(std::move(s)), a(s.a) {}

    Scalar<MultiPrecision, true>::Scalar(const Scalar<MultiPrecision, false>& s)
            : Scalar<MultiPrecision, false>(s), a(0) {}

    Scalar<MultiPrecision, true>::Scalar(Scalar<MultiPrecision, false>&& s) noexcept
            : Scalar<MultiPrecision, false>(std::move(s)), a(0) {}

    Scalar<MultiPrecision, true>::Scalar(int i, ScalarUnit a)
            : Scalar<MultiPrecision, false>(i), a(a) {}

    Scalar<MultiPrecision, true>::Scalar(SignedScalarUnit unit, ScalarUnit a)
            : Scalar<MultiPrecision, false>(unit), a(a) {}

    Scalar<MultiPrecision, true>::Scalar(double d, ScalarUnit a)
            : Scalar<MultiPrecision, false>(d), a(a) {}

    Scalar<MultiPrecision, true>::Scalar(const char* s, ScalarUnit a)
            : Scalar<MultiPrecision, false>(s), a(a) {}

    Scalar<MultiPrecision, true>::Scalar(const wchar_t* s, ScalarUnit a)
            : Scalar<MultiPrecision, false>(s), a(a) {}

    Scalar<MultiPrecision, true>& Scalar<MultiPrecision, true>::operator=(
            const Scalar<MultiPrecision, true>& s) {
        a = s.a;
        Scalar<MultiPrecision, false>::operator=(s);
    }

    Scalar<MultiPrecision, true>& Scalar<MultiPrecision, true>::operator=(
            Scalar<MultiPrecision, true>&& s) noexcept {
        a = s.a;
        Scalar<MultiPrecision, false>::operator=(std::move(s));
    }

    Scalar<MultiPrecision, true>& Scalar<MultiPrecision, true>::operator=(const Scalar<MultiPrecision, false>& s) {
        a = 0;
        Scalar<MultiPrecision, false>::operator=(s);
        return *this;
    }

    Scalar<MultiPrecision, true>& Scalar<MultiPrecision, true>::operator=(Scalar<MultiPrecision, false>&& s) noexcept {
        a = 0;
        Scalar<MultiPrecision, false>::operator=(std::move(s));
        return *this;
    }

    Scalar<MultiPrecision, true> Scalar<MultiPrecision, true>::operator+(
            const Scalar<MultiPrecision, false>& s) const {
        auto result = add<true>(*this, s);
        cutLength(result);
        if(getA() != 0)
            result.applyError(getAccuracy());
        return result;
    }

    Scalar<MultiPrecision, true> Scalar<MultiPrecision, true>::operator-(
            const Scalar<MultiPrecision, false>& s) const {
        Q_UNUSED(MultiPrecision)
        auto result = sub<true>(*this, s);
        cutLength(result);
        if(getA() != 0)
            result.applyError(getAccuracy());
        return result;
    }

    Scalar<MultiPrecision, true> Scalar<MultiPrecision, true>::operator*(
            const Scalar<MultiPrecision, false>& s) const {
        Q_UNUSED(MultiPrecision)
        auto result = mul<true>(*this, s);
        cutLength(result);
        if(getA() != 0)
            result.applyError(s * getAccuracy());
        return result;
    }

    Scalar<MultiPrecision, true> Scalar<MultiPrecision, true>::operator/(
            const Scalar<MultiPrecision, false>& s) const {
        Q_UNUSED(MultiPrecision)
        auto result = div<true>(*this, s);
        if(getA() != 0)
            result.applyError(getAccuracy() / s);
        return result;
    }

    Scalar<MultiPrecision, true> Scalar<MultiPrecision, true>::operator+(
            const Scalar<MultiPrecision, true>& s) const {
        auto result = add<true>(*this, s);
        cutLength(result);
        if(getA() != 0 || s.getA() != 0) {
            if(getA() == 0)
                result.applyError(s.getAccuracy());
            else if(s.getA() == 0)
                result.applyError(getAccuracy());
            else
                result.applyError(getAccuracy() + s.getAccuracy());
        }
        return result;
    }

    Scalar<MultiPrecision, true> Scalar<MultiPrecision, true>::operator-(
            const Scalar<MultiPrecision, true>& s) const {
        auto result = sub<true>(*this, s);
        cutLength(result);
        if(getA() != 0 || s.getA() != 0) {
            if(getA() == 0)
                result.applyError(s.getAccuracy());
            else if(s.getA() == 0)
                result.applyError(getAccuracy());
            else
                result.applyError(getAccuracy() + s.getAccuracy());
        }
        return result;
    }

    Scalar<MultiPrecision, true> Scalar<MultiPrecision, true>::operator*(
            const Scalar<MultiPrecision, true>& s) const {
        auto result = mul<true>(*this, s);
        cutLength(result);
        if(getA() != 0 || s.getA() != 0) {
            const Scalar<MultiPrecision, false>& casted = s;
            if(getA() == 0)
                result.applyError(static_cast<const Scalar<MultiPrecision, false>&>(*this) * s.getAccuracy());
            else if(s.getA() == 0)
                result.applyError(casted * getAccuracy());
            else {
                Scalar<MultiPrecision, false> this_a = getAccuracy();
                Scalar<MultiPrecision, false> s_a = s.getAccuracy();
                Scalar<MultiPrecision, false> temp1 = this_a * s_a;
                Scalar<MultiPrecision, false> temp2
                        = static_cast<const Scalar<MultiPrecision, false>&>(*this) * s_a + casted * this_a;
                result.applyError(temp1 + temp2);
            }
        }
        return result;
    }

    Scalar<MultiPrecision, true> Scalar<MultiPrecision, true>::operator/(
            const Scalar<MultiPrecision, true>& s) const {
        auto result = div<true>(*this, s);
        if(getA() != 0 || s.getA() != 0) {
            const Scalar<MultiPrecision, false>& casted = s;
            if(getA() == 0) {
                Scalar<MultiPrecision, false> s_a = s.getAccuracy();
                Scalar<MultiPrecision, false> temp_1
                        = static_cast<const Scalar<MultiPrecision, false>&>(*this) * s_a;
                Scalar<MultiPrecision, false> temp_2 = casted * (casted - s_a);
                result.applyError(temp_1 / temp_2);
            }
            else if(s.getA() == 0)
                result.applyError(getAccuracy() / s);
            else {
                Scalar<MultiPrecision, false> s_a = s.getAccuracy();
                Scalar<MultiPrecision, false> temp_1
                        = static_cast<const Scalar<MultiPrecision, false>&>(*this) * s_a + casted * getAccuracy();
                Scalar<MultiPrecision, false> temp_2 = casted * (casted - s_a);
                result.applyError(temp_1 / temp_2);
            }
        }
        return result;
    }

    Scalar<MultiPrecision, true> Scalar<MultiPrecision, true>::operator<<(int bits) const {
        if(bits == 0)
            return Scalar(*this);
        if(bits < 0)
            return *this >> -bits;
        int size = Scalar<MultiPrecision, false>::getSize();
        const int quotient = bits / ScalarUnitWidth; //NOLINT: quotient < INT_MAX
        const unsigned int remainder = bits - quotient * ScalarUnitWidth;
        //If remainder = 0, we must return directly because shifting a ScalarUnit for ScalarUnitWidth bits is a undefined behavior.
        if(remainder == 0) {
            Scalar copy(*this);
            copy.power += quotient;
            return copy;
        }

        const int size_1 = size - 1;
        const bool carry = countLeadingZeros(byte[size_1]) < remainder;
        const bool a_carry = countLeadingZeros(a) < remainder;
        Scalar result(length >= 0 ? (size + carry - a_carry) : -(size + carry - a_carry), power + quotient + carry,
                      a_carry ? a >> (ScalarUnitWidth - remainder) : a << remainder); //Add 1 to a if a_carry is true to get more accurate estimation.
        const auto byte_start = byte + a_carry;
        size -= a_carry;
        result[0] = 0;
        for(int i = 0; i < size_1; ++i) {
            result[i] |= byte_start[i] << remainder;
            result[i + 1] = byte_start[i] >> (ScalarUnitWidth - remainder);
        }
        result[size_1] |= byte_start[size_1] << remainder;
        if(carry)
            result[size] = byte_start[size_1] >> (ScalarUnitWidth - remainder);
        return result;
    }

    Scalar<MultiPrecision, true> Scalar<MultiPrecision, true>::operator>>(int bits) const {
        if(bits == 0)
            return Scalar(*this);
        if(bits < 0)
            return *this << -bits;
        const int size = Scalar<MultiPrecision, false>::getSize();
        const int quotient = bits / ScalarUnitWidth; //NOLINT: quotient < INT_MAX
        const unsigned int remainder = bits - quotient * ScalarUnitWidth;
        //If remainder = 0, we must return directly because shifting a ScalarUnit for ScalarUnitWidth bits is a undefined behavior.
        if(remainder == 0) {
            Scalar copy(*this);
            copy.power -= quotient;
            return copy;
        }

        const int size_1 = size - 1;
        const bool carry = (countLeadingZeros(byte[size_1]) + remainder) < ScalarUnitWidth;
        const ScalarUnit accuracy = a >> remainder;
        Scalar result(length >= 0 ? (size + carry) : -(size + carry)
                , power - quotient + carry - 1, accuracy > 0 ? accuracy : 0);
        if(carry)
            result[size] = byte[size_1] >> remainder;

        for(int i = size_1; i > 0; --i) {
            result[i] = byte[i] << (ScalarUnitWidth - remainder);
            result[i] |= byte[i - 1] >> remainder;
        }
        result[0] = byte[0] << (ScalarUnitWidth - remainder);
        return result;
    }

    Scalar<MultiPrecision, true> Scalar<MultiPrecision, true>::operator-() const {
        Scalar result(-length, power, a);
        memcpy(result.byte, byte, Scalar<MultiPrecision, false>::getSize() * sizeof(ScalarUnit));
        return result;
    }
    /*!
     * Add error to this and adjust this->length as well as this->byte.
     *
     * Optimize:
     * error is always non-zero, if(!error.isZero()) is unnecessary.
     */
    Scalar<MultiPrecision, true>& Scalar<MultiPrecision, true>::applyError(
            const Scalar<MultiPrecision, false>& error) {
        if(!error.isZero()) {
            int size = Scalar<MultiPrecision, false>::getSize();
            int copy = size;
            int temp = power - size + 1 - error.power;
            ScalarUnit copy_a = a;
            if(temp <= 0) {
                if(temp < 0) {
                    auto error_1 = getAccuracy() + error;
                    size += temp;
                    //Use (a += error_1.byte[error_1.getSize() - 1] + 1) for more conservative error estimate.
                    a += error_1.byte[error_1.getSize() - 1];
                }
                else
                    //Use (a += error.byte[error.getSize() - 1] + 1) for more conservative error estimate.
                    a += error.byte[error.getSize() - 1];
            }

            if(a < copy_a) {
                //Use (a = 2) for more conservative error estimate.
                a = 1;
                --size;
            }

            if(size < 1) {
                qWarning("Accumulated too many errors.");
                power += 1 - size;
                length = length < 0 ? -1 : 1;
                return *this;
            }

            if(size < copy) {
                auto new_byte = reinterpret_cast<ScalarUnit*>(malloc(size * sizeof(ScalarUnit)));
                memcpy(new_byte, byte + copy - size, size * sizeof(ScalarUnit));
                this->~Scalar();
                byte = new_byte;
            }
            length = length < 0 ? -size : size;
        }
        return *this;
    }

    void Scalar<MultiPrecision, true>::swap(
            Scalar<MultiPrecision, true>& s) noexcept {
        std::swap(a, s.a);
        Scalar<MultiPrecision, false>::swap(s);
    }
    //!Return accuracy in class Scalar.
    Scalar<MultiPrecision, false> Scalar<MultiPrecision, true>::getAccuracy() const {
        Scalar<MultiPrecision, false> result(
                1, power - Scalar<MultiPrecision, true>::getSize() + 1);
        result[0] = a;
        return result;
    }
    ///////////////////////////////////////////Float-Double////////////////////////////////////////////////
    /*!
     * Inspect if a float is a integer through its binary expression.
     */
    bool Scalar<Float, false>::isInteger() const {
        float_extract extract{f};
        const auto exp = extract.exp;
        if(exp == 0U)
            return extract.value == 0;

        unsigned int zeros;
        if(extract.low == 0U) {
            if(extract.high == 0U)
                return true;
            else
                zeros = countBackZeros(extract.high) + 16; //extract.low is zero, which has 16 bits.
        }
        else
            zeros = countBackZeros(extract.low);
        /*
         * exp + zeros - 127 >= 23
         * , 127 is the exp bias of float numbers, 23 is the number of bits of significand of float numbers.
         * We move 127 to the right hand side to avoid underflow.
         */
        return exp + zeros >= 150;
    }

    void Scalar<Float, true>::swap(Scalar<Float, true>& s) noexcept {
        std::swap(a, s.a);
        Scalar<Float, false>::swap(s);
    }
    /*!
     * Inspect if a float is a integer through its binary expression.
     */
    bool Scalar<Double, false>::isInteger() const {
        double_extract extract{d};
        const auto exp = extract.exp;
        if(exp == 0U)
            return extract.value == 0;

        unsigned int zeros;
        if(extract.low == 0U) {
            if(extract.high == 0U)
                return true;
            else
                zeros = countBackZeros(extract.high) + 32; //extract.low is zero, which has 32 bits.
        }
        else
            zeros = countBackZeros(extract.low);
        /*
         * exp + zeros - 1023 >= 52
         * , 1023 is the exp bias of float numbers, 52 is the number of bits of significand of float numbers.
         * We move 1023 to the right hand side to avoid underflow.
         */
        return exp + zeros >= 1075;
    }

    void Scalar<Double, true>::swap(Scalar<Double, true>& s) noexcept {
        std::swap(a, s.a);
        Scalar<Double, false>::swap(s);
    }
}