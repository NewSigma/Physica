/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_MULTISCALAR_H
#define PHYSICA_MULTISCALAR_H

#ifndef PHYSICA_SCALAR_H
    #include "Core/Header/MultiPrecition/Scalar.h"
#endif

#include <iomanip>
#include <Core/Header/MultiPrecition/Scalar.h>

#include "Core/MultiPrecision/Basic/AddBasic.h"
#include "Core/MultiPrecision/Basic/DivBasic.h"
#include "Core/MultiPrecision/Util/ArraySupport.h"
#include "Core/Header/Utils/MetaSupport.h"

namespace Physica::Core {
    ///////////////////////////////////////////MultiPrecision/////////////////////////////////////////
    //!Reference: MaTHmu Project Group.计算机代数系统的数学原理[M].Beijing: TsingHua University Press, 2009.45
    template<size_t maxPrecision, bool errorTrack>
    inline Scalar<maxPrecision, errorTrack> operator^(
            const Scalar<maxPrecision, errorTrack>& s1,
            const Scalar<maxPrecision, errorTrack>& s2) {
        return s1.isInteger() ? powScalar(s1, s2) : exp(ln(s1) * s2);
    }

    template<size_t maxPrecision, bool errorTrack>
    inline Scalar<maxPrecision, errorTrack>& operator++(
            Scalar<maxPrecision, errorTrack>& s) {
        s += BasicConst::getInstance().get_1();
        return s;
    }

    template<size_t maxPrecision, bool errorTrack>
    inline Scalar<maxPrecision, errorTrack>& operator--(
            Scalar<maxPrecision, errorTrack>& s) {
        s -= BasicConst::getInstance().get_1();
        return s;
    }

    template<size_t maxPrecision, bool errorTrack>
    inline Scalar<maxPrecision, errorTrack> operator++( //NOLINT confusing-warning
            Scalar<maxPrecision, errorTrack>& s, int) {
        Scalar<maxPrecision, errorTrack> temp(s);
        s += BasicConst::getInstance().get_1();
        return temp;
    }

    template<size_t maxPrecision, bool errorTrack>
    inline Scalar<maxPrecision, errorTrack> operator--( //NOLINT confusing-warning
            Scalar<maxPrecision, errorTrack>& s, int) {
        Scalar<maxPrecision, errorTrack> temp(s);
        s -= BasicConst::getInstance().get_1();
        return temp;
    }
    //////////////////////////////////MultiPrecision-WithoutError///////////////////////////////////
    template<size_t maxPrecision>
    Scalar<maxPrecision, false>::Scalar() noexcept
            : byte(nullptr), length(0), power(0) { Q_UNUSED(maxPrecision) }

    template<size_t maxPrecision>
    Scalar<maxPrecision, false>::Scalar(int length, int power) noexcept
            : byte(reinterpret_cast<ScalarUnit*>(malloc(abs(length) * sizeof(ScalarUnit))))
            , length(length), power(power) { Q_UNUSED(maxPrecision) }

    template<size_t maxPrecision>
    Scalar<maxPrecision, false>::Scalar(const Scalar<maxPrecision, false>& s)
            : byte(reinterpret_cast<ScalarUnit*>(malloc(s.getSize() * sizeof(ScalarUnit))))
            , length(s.length), power(s.power) {
        memcpy(byte, s.byte, getSize() * sizeof(ScalarUnit));
    }

    template<size_t maxPrecision>
    Scalar<maxPrecision, false>::Scalar(Scalar<maxPrecision, false>&& s) noexcept
            : byte(s.byte), length(s.length), power(s.power) {
        s.byte = nullptr;
    }

    template<size_t maxPrecision>
    Scalar<maxPrecision, false>::Scalar(SignedScalarUnit unit)
            : byte(reinterpret_cast<ScalarUnit*>(malloc(sizeof(ScalarUnit))))
            , length(unit > 0 ? length : -length), power(0) {
        Q_UNUSED(maxPrecision)
        byte[0] = unit > 0 ? unit : -unit;
    }

    template<size_t maxPrecision>
    Scalar<maxPrecision, false>::Scalar(double d) {
        Q_UNUSED(maxPrecision)
        if(d == 0) {
            byte = reinterpret_cast<ScalarUnit*>(malloc(sizeof(ScalarUnit)));
            length = 1;
            byte[0] = power = 0;
            return;
        }
        double_extract extract{d};
        auto quotient = static_cast<int>(extract.structure.exp) - 1023;
        power = quotient / __WORDSIZE;
        //Have power * __WORDSIZE < extract.structure.exp, so that remainder > 0.
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
                byte[1] += static_cast<ScalarUnit>(extract.structure.high) >> (20 - remainder);
                byte[0] = static_cast<ScalarUnit>(extract.structure.high) << (44 + remainder);
                byte[0] += static_cast<ScalarUnit>(extract.structure.low) << (32 - (20 - remainder));
            }
            else {
                byte[1] += static_cast<ScalarUnit>(extract.structure.high) << (remainder - 20);
                byte[1] += static_cast<ScalarUnit>(extract.structure.low) >> (32 - (remainder - 20));
                byte[0] = static_cast<ScalarUnit>(extract.structure.low) << (32 + (remainder - 20));
            }
        }
        else {
            length = 1;
            byte = reinterpret_cast<ScalarUnit*>(malloc(sizeof(ScalarUnit)));
            //Hidden bit
            byte[0] = 1;
            byte[0] <<= 20U;
            byte[0] += static_cast<ScalarUnit>(extract.structure.high);
            byte[0] <<= 32U;
            byte[0] += static_cast<ScalarUnit>(extract.structure.low);
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
            byte[2] += static_cast<ScalarUnit>(extract.structure.high) >> (20 - remainder);
            byte[1] = static_cast<ScalarUnit>(extract.structure.high) << (32 - (20 - remainder));
            byte[1] +=  static_cast<ScalarUnit>(extract.structure.low) >> (20 - remainder);
            byte[0] = static_cast<ScalarUnit>(extract.structure.low) << remainder;
        }
        else {
            length = 2;
            byte = reinterpret_cast<ScalarUnit*>(malloc(length * sizeof(ScalarUnit)));
            //Hidden bit
            byte[1] = 1;
            byte[1] <<= remainder;
            byte[1] += static_cast<ScalarUnit>(extract.structure.high) << (remainder - 20);
            byte[1] += static_cast<ScalarUnit>(extract.structure.low) >> (32 - (remainder - 20));
            byte[0] = static_cast<ScalarUnit>(extract.structure.low) << (remainder - 20);
        }
    #endif
        if(extract.structure.sign)
            length = -length;
    }
    /*!
     * Not accurate.
     */
    template<size_t maxPrecision>
    Scalar<maxPrecision, false>::Scalar(const char* s)
            : Scalar(strtod(s, nullptr)) { Q_UNUSED(maxPrecision) }
    /*!
     * Not accurate.
     */
    template<size_t maxPrecision>
    Scalar<maxPrecision, false>::Scalar(const wchar_t* s) {
        Q_UNUSED(maxPrecision)
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

    template<size_t maxPrecision>
    Scalar<maxPrecision, false>::~Scalar() {  Q_UNUSED(maxPrecision) free(byte); }

    template<size_t maxPrecision>
    Scalar<maxPrecision, false>& Scalar<maxPrecision, false>::operator= (
            const Scalar<maxPrecision, false>& s) {
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

    template<size_t maxPrecision>
    Scalar<maxPrecision, false>& Scalar<maxPrecision, false>::operator=(
            Scalar<maxPrecision, false>&& s) noexcept {
        this->~Scalar();
        byte = s.byte;
        s.byte = nullptr;
        length = s.length;
        power = s.power;
        return *this;
    }

    template<size_t maxPrecision>
    Scalar<maxPrecision, false>::operator double() const {
        Q_UNUSED(maxPrecision)
        if(isZero())
            return 0.0;
        double_extract extract{0};
        extract.structure.sign = length < 0;

        const auto zeroCount = countLeadingZeros(byte[getSize() - 1]);
        auto exp = power * __WORDSIZE + ScalarUnitWidth - zeroCount - 1 + 1023;
        if(exp >= 2047) {
            extract.structure.high = extract.structure.low = 0;
            extract.structure.exp = 2047;
            return extract.value;
        }
        if(exp <= 0) {
            return 0.0;
        }
        extract.structure.exp = exp;

        auto size = getSize();
        auto temp = byte[size - 1] << (zeroCount + 1);
    #ifdef PHYSICA_64BIT
        extract.structure.high = temp >> 44U;
        if(zeroCount <= 11) {
            extract.structure.low = temp << 20U >> 32U;
        }
        else {
            if(zeroCount <= 44 - 1) {
                extract.structure.low = temp << 20U >> 32U;
                if(size > 1)
                    extract.structure.low += byte[size - 2] >> (64 - (32 - (64 - 20 - zeroCount - 1)));
            }
            else {
                if(size > 1) {
                    extract.structure.high += byte[size - 2] >> (64 - (20 - (64 - zeroCount - 1)));
                    extract.structure.low = byte[size - 2] << (20 - (64 - zeroCount - 1)) >> 32U;
                }
            }
        }
    #endif
    #ifdef PHYSICA_32BIT
        extract.structure.high = temp >> 12U;
        if(zeroCount <= 11) {
            extract.structure.low = temp << 20U;
            if(size > 1)
                extract.structure.low = byte[size - 1] >> (32 - 20 - zeroCount - 1);
        }
        else {
            if(size > 1) {
                extract.structure.high += byte[size - 1] >> (32 - (zeroCount + 1 - 12));
                extract.structure.low = byte[size - 1] << (zeroCount + 1 - 12);
            }
            if(size > 2)
                extract.structure.low += byte[size - 2] >> (32 - (zeroCount + 1 - 12));
        }
    #endif
        return extract.value;
    }

    template<size_t maxPrecision>
    template<size_t maxPrecision2>
    Scalar<META_MAX<maxPrecision, maxPrecision2>::value, false> Scalar<maxPrecision, false>::operator+(
            const Scalar<maxPrecision2, false>& s) const {
        auto result = add(*this, s);
        cutLength(result);
        return result;
    }

    template<size_t maxPrecision>
    template<size_t maxPrecision2>
    Scalar<META_MAX<maxPrecision, maxPrecision2>::value, false> Scalar<maxPrecision, false>::operator-(
            const Scalar<maxPrecision2, false>& s) const {
        auto result = sub(*this, s);
        cutLength(result);
        return result;
    }

    template<size_t maxPrecision>
    template<size_t maxPrecision2>
    Scalar<META_MAX<maxPrecision, maxPrecision2>::value, false> Scalar<maxPrecision, false>::operator*(
            const Scalar<maxPrecision2, false>& s) const {
        auto result = mul(*this, s);
        cutLength(result);
        return result;
    }

    template<size_t maxPrecision>
    template<size_t maxPrecision2>
    Scalar<META_MAX<maxPrecision, maxPrecision2>::value, false> Scalar<maxPrecision, false>::operator/(
            const Scalar<maxPrecision2, false>& s) const {
        return div(*this, s);
    }

    template<size_t maxPrecision>
    Scalar<maxPrecision, false> Scalar<maxPrecision, false>::operator<<(int bits) const {
        Q_UNUSED(maxPrecision)
        if(bits == 0)
            return Scalar(*this);
        if(bits < 0)
            return *this >> -bits;
        const int size = getSize();
        const int quotient = bits / ScalarUnitWidth; //NOLINT: quotient < INT_MAX
        const unsigned int remainder = bits - quotient * ScalarUnitWidth;
        const bool carry = countLeadingZeros(byte[size - 1]) < remainder;

        Scalar result(length >= 0 ? (size + carry) : -(size + carry), power + quotient + carry);
        result[0] = 0;
        for(int i = 0; i < size - 1; ++i) {
            result[i] |= byte[i] << remainder;
            result[i + 1] = byte[i] >> (ScalarUnitWidth - remainder);
        }
        result[size - 1] |= byte[size - 1] << remainder;
        if(carry)
            result[size] = byte[size - 1] >> (ScalarUnitWidth - remainder);
        return result;
    }

    template<size_t maxPrecision>
    Scalar<maxPrecision, false> Scalar<maxPrecision, false>::operator>>(int bits) const {
        Q_UNUSED(maxPrecision)
        if(bits == 0)
            return Scalar(*this);
        if(bits < 0)
            return *this << -bits;
        const int size = getSize();
        const int quotient = bits / ScalarUnitWidth; //NOLINT: quotient < INT_MAX
        const unsigned int remainder = bits - quotient * ScalarUnitWidth;
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

    template<size_t maxPrecision>
    Scalar<maxPrecision, false> Scalar<maxPrecision, false>::operator-() const {
        Q_UNUSED(maxPrecision)
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
    template<size_t maxPrecision>
    bool absCompare(const Scalar<maxPrecision, false>& s1
            , const Scalar<maxPrecision, false>& s2) {
        if(s1.getPower() > s2.getPower())
            return true;
        if(s1.getPower() < s2.getPower())
            return false;
        const Scalar<maxPrecision, false>* longer, *shorter;
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
    template<size_t maxPrecision>
    bool operator> (const Scalar<maxPrecision, false>& s1
            , const Scalar<maxPrecision, false>& s2) {
        //Judge from sign.
        if(s1.isPositive()) {
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
            result = true;
        else if(s1.getPower() < s2.getPower())
            result = false;
        else {
            //The only method left.
            Scalar<maxPrecision, false> subtract = s1 - s2;
            result = subtract.isPositive();
        }
        return result;
    }
    /*!
     * Optimize:
     * Is subtract faster than comparing the elements?
     */
    template<size_t maxPrecision>
    bool operator< (const Scalar<maxPrecision, false>& s1
            , const Scalar<maxPrecision, false>& s2) {
        //Judge from sign.
        if(s1.isPositive()) {
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
            result = false;
        else if(s1.getPower() < s2.getPower())
            result = true;
        else {
            //The only method left.
            Scalar<maxPrecision, false> subtract = s1 - s2;
            result = subtract.isNegative();
        }
        return result;
    }

    template<size_t maxPrecision>
    bool operator>= (const Scalar<maxPrecision, false>& s1
            , const Scalar<maxPrecision, false>& s2) {
        return !(s1 < s2);
    }

    template<size_t maxPrecision>
    bool operator<= (const Scalar<maxPrecision, false>& s1
            , const Scalar<maxPrecision, false>& s2) {
        return !(s1 > s2);
    }

    template<size_t maxPrecision>
    bool operator== (const Scalar<maxPrecision, false>& s1
            , const Scalar<maxPrecision, false>& s2) {
        return s1.getPower() == s2.getPower()
               //Here length may not equal n.length because we define numbers like 1.0 and 1.00 are equal.
               && ((s1.getLength() ^ s2.getLength()) >= 0) //NOLINT
               && s1.getA() == s2.getA()
               && (s1 - s2).isZero();
    }

    template<size_t maxPrecision>
    bool operator!= (const Scalar<maxPrecision, false>& s1
            , const Scalar<maxPrecision, false>& s2) {
        return !(s1 == s2);
    }

    template<size_t maxPrecision>
    void Scalar<maxPrecision, false>::swap(
            Scalar<maxPrecision, false>& s) noexcept {
        std::swap(byte, s.byte);
        std::swap(length, s.length);
        std::swap(power, s.power);
    }
    ///////////////////////////////////MultiPrecision-WithError/////////////////////////////////////
    template<size_t maxPrecision>
    Scalar<maxPrecision, true>::Scalar() noexcept
            : Scalar<maxPrecision, false>(), a(0) { Q_UNUSED(maxPrecision) }

    template<size_t maxPrecision>
    Scalar<maxPrecision, true>::Scalar(int length, int power, ScalarUnit a) noexcept
            : Scalar<maxPrecision, false>(length, power), a(a) { Q_UNUSED(maxPrecision) }

    template<size_t maxPrecision>
    Scalar<maxPrecision, true>::Scalar(const Scalar<maxPrecision, true>& s)
            : Scalar<maxPrecision, false>(s), a(s.a) {}

    template<size_t maxPrecision>
    Scalar<maxPrecision, true>::Scalar(Scalar<maxPrecision, true>&& s) noexcept
            : Scalar<maxPrecision, false>(std::move(s))
            , a(s.a) {}

    template<size_t maxPrecision>
    Scalar<maxPrecision, true>::Scalar(SignedScalarUnit unit, ScalarUnit a)
            : Scalar<maxPrecision, false>(unit), a(a) { Q_UNUSED(maxPrecision) }

    template<size_t maxPrecision>
    Scalar<maxPrecision, true>::Scalar(double d, ScalarUnit a)
            : Scalar<maxPrecision, false>(d), a(a) { Q_UNUSED(maxPrecision) }

    template<size_t maxPrecision>
    Scalar<maxPrecision, true>::Scalar(const char* s, ScalarUnit a)
            : Scalar<maxPrecision, false>(s), a(a) { Q_UNUSED(maxPrecision) }

    template<size_t maxPrecision>
    Scalar<maxPrecision, true>::Scalar(const wchar_t* s, ScalarUnit a)
            : Scalar<maxPrecision, false>(s), a(a) { Q_UNUSED(maxPrecision) }

    template<size_t maxPrecision>
    Scalar<maxPrecision, true>& Scalar<maxPrecision, true>::operator=(
            const Scalar<maxPrecision, true>& s) {
        a = s.a;
        Scalar<maxPrecision, false>::operator=(s);
    }

    template<size_t maxPrecision>
    Scalar<maxPrecision, true>& Scalar<maxPrecision, true>::operator=(
            Scalar<maxPrecision, true>&& s) noexcept {
        a = s.a;
        Scalar<maxPrecision, false>::operator=(std::move(s));
    }

    template<size_t maxPrecision>
    template<size_t maxPrecision2>
    Scalar<META_MAX<maxPrecision, maxPrecision2>::value, true> Scalar<maxPrecision, true>::operator+(
            const Scalar<maxPrecision2, true>& s) const {
        auto result = Scalar<maxPrecision, false>::operator+(s);
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

    template<size_t maxPrecision>
    template<size_t maxPrecision2>
    Scalar<META_MAX<maxPrecision, maxPrecision2>::value, true> Scalar<maxPrecision, true>::operator-(
            const Scalar<maxPrecision2, true>& s) const {
        auto result = Scalar<maxPrecision, false>::operator-(s);
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

    template<size_t maxPrecision>
    template<size_t maxPrecision2>
    Scalar<META_MAX<maxPrecision, maxPrecision2>::value, true> Scalar<maxPrecision, true>::operator*(
            const Scalar<maxPrecision2, true>& s) const {
        auto result = Scalar<maxPrecision, false>::operator*(s);
        if(getA() != 0 || s.getA() != 0) {
            if(getA() == 0)
                result.applyError(*this * s.getAccuracy());
            else if(s.getA() == 0)
                result.applyError(s * getAccuracy());
            else {
                Scalar this_a = getAccuracy();
                Scalar s_a = s.getAccuracy();
                Scalar temp1 = this_a * s_a;
                Scalar temp2 = *this * s_a + s * this_a;
                result.applyError(temp1 + temp2);
            }
        }
        return result;
    }

    template<size_t maxPrecision>
    template<size_t maxPrecision2>
    Scalar<META_MAX<maxPrecision, maxPrecision2>::value, true> Scalar<maxPrecision, true>::operator/(
            const Scalar<maxPrecision2, true>& s) const {
        auto result = Scalar<maxPrecision, false>::operator/(s);
        if(getA() != 0 || s.getA() != 0) {
            if(getA() == 0) {
                Scalar s_a = s.getAccuracy();
                Scalar temp_1 = *this * s_a;
                Scalar temp_2 = s * (s - s_a);
                result.applyError(temp_1 / temp_2);
            }
            else if(s.getA() == 0)
                result.applyError(getAccuracy() / s);
            else {
                Scalar s_a = s.getAccuracy();
                Scalar temp_1 = *this * s_a + s * getAccuracy();
                Scalar temp_2 = s * (s - s_a);
                result.applyError(temp_1 / temp_2);
            }
        }
        return result;
    }

    template<size_t maxPrecision>
    Scalar<maxPrecision, true> Scalar<maxPrecision, true>::operator<<(int bits) const {
        if(bits == 0)
            return Scalar(*this);
        if(bits < 0)
            return *this >> -bits;
        int size = Scalar<maxPrecision, false>::getSize();
        const int quotient = bits / ScalarUnitWidth; //NOLINT: quotient < INT_MAX
        const unsigned int remainder = bits - quotient * ScalarUnitWidth;
        const bool carry = countLeadingZeros(byte[size - 1]) < remainder;
        const bool a_carry = countLeadingZeros(a) < remainder;

        Scalar result(length >= 0 ? (size + carry - a_carry) : -(size + carry - a_carry), power + quotient + carry,
                      a_carry ? a >> (ScalarUnitWidth - remainder) : a << remainder); //Add 1 to a if a_carry is true to get more accurate estimation.
        const auto byte_start = byte + a_carry;
        size -= a_carry;
        result[0] = 0;
        for(int i = 0; i < size - 1; ++i) {
            result[i] |= byte_start[i] << remainder;
            result[i + 1] = byte_start[i] >> (ScalarUnitWidth - remainder);
        }
        result[size - 1] |= byte_start[size - 1] << remainder;
        if(carry)
            result[size] = byte_start[size - 1] >> (ScalarUnitWidth - remainder);
        return result;
    }

    template<size_t maxPrecision>
    Scalar<maxPrecision, true> Scalar<maxPrecision, true>::operator>>(int bits) const {
        if(bits == 0)
            return Scalar(*this);
        if(bits < 0)
            return *this << -bits;
        const int size = Scalar<maxPrecision, false>::getSize();
        const int quotient = bits / ScalarUnitWidth; //NOLINT: quotient < INT_MAX
        const unsigned int remainder = bits - quotient * ScalarUnitWidth;
        const bool carry = (countLeadingZeros(byte[size - 1]) + remainder) < ScalarUnitWidth;
        const ScalarUnit accuracy = a >> remainder;

        Scalar result(length >= 0 ? (size + carry) : -(size + carry)
                , power - quotient + carry - 1, accuracy > 0 ? accuracy : 0);
        if(carry)
            result[size] = byte[size - 1] >> remainder;

        for(int i = size - 1; i > 0; --i) {
            result[i] = byte[i] << (ScalarUnitWidth - remainder);
            result[i] |= byte[i - 1] >> remainder;
        }
        result[0] = byte[0] << (ScalarUnitWidth - remainder);
        return result;
    }

    template<size_t maxPrecision>
    Scalar<maxPrecision, true> Scalar<maxPrecision, true>::operator-() const {
        Scalar result(-length, power, a);
        memcpy(result.byte, byte, Scalar<maxPrecision, false>::getSize() * sizeof(ScalarUnit));
        return result;
    }
    /*!
     * Add error to this and adjust this->length as well as this->byte.
     *
     * Optimize:
     * error is always non-zero, if(!error.isZero()) is unnecessary.
     */
    template<size_t maxPrecision>
    Scalar<maxPrecision, true>& Scalar<maxPrecision, true>::applyError(
            const Scalar<maxPrecision, true>& error) {
        if(!error.isZero()) {
            int size = Scalar<maxPrecision, false>::getSize();
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

    template<size_t maxPrecision>
    inline void Scalar<maxPrecision, true>::swap(
            Scalar<maxPrecision, true>& s) noexcept {
        std::swap(a, s.a);
        Scalar<maxPrecision, false>::swap(s);
    }
    //!Return accuracy in class Scalar.
    template<size_t maxPrecision>
    Scalar<maxPrecision, false> Scalar<maxPrecision, true>::getAccuracy() const {
        Scalar<maxPrecision, false> result(
                1, power - Scalar<maxPrecision, true>::getSize() + 1);
        result[0] = a;
        return result;
    }

    template<size_t maxPrecision>
    inline Scalar<maxPrecision, true> Scalar<maxPrecision, true>::getMaximum() const {
        Q_UNUSED(maxPrecision)
        return *this + getAccuracy();
    }

    template<size_t maxPrecision>
    inline Scalar<maxPrecision, true> Scalar<maxPrecision, true>::getMinimum() const {
        Q_UNUSED(maxPrecision)
        return *this - getAccuracy();
    }
    ///////////////////////////////////////////Float-Double////////////////////////////////////////////////
    /////////////////////////////////////////////Float////////////////////////////////////////////////
    template<bool errorTrack>
    inline Scalar<0, errorTrack>& operator++(
            Scalar<0, errorTrack>& s) {
        s += 1.0F;
        return s;
    }

    template<bool errorTrack>
    inline Scalar<0, errorTrack>& operator--(
            Scalar<0, errorTrack>& s) {
        s -= 1.0F;
        return s;
    }

    template<bool errorTrack>
    inline Scalar<0, errorTrack> operator++( //NOLINT confusing-warning
            Scalar<0, errorTrack>& s, int) {
        Scalar<0, errorTrack> temp(s);
        s += 1.0F;
        return temp;
    }

    template<bool errorTrack>
    inline Scalar<0, errorTrack> operator--( //NOLINT confusing-warning
            Scalar<0, errorTrack>& s, int) {
        Scalar<0, errorTrack> temp(s);
        s -= 1.0F;
        return temp;
    }
    ////////////////////////////////////////Float-WithoutError///////////////////////////////////////////
    inline Scalar<0, false>::Scalar() : f(0) {}

    inline Scalar<0, false>::Scalar(float f) : f(f) {}
    /////////////////////////////////////////Float-WithError////////////////////////////////////////////////
    inline Scalar<0, true>::Scalar() : Scalar<0, false>(), a(0) {}

    inline Scalar<0, true>::Scalar(float f, float a) : Scalar<0, false>(f), a(a) {}
    /////////////////////////////////////////////Double////////////////////////////////////////////////
    template<bool errorTrack>
    inline Scalar<1, errorTrack>& operator++(
            Scalar<1, errorTrack>& s) {
        s += 1.0;
        return s;
    }

    template<bool errorTrack>
    inline Scalar<1, errorTrack>& operator--(
            Scalar<1, errorTrack>& s) {
        s -= 1.0;
        return s;
    }

    template<bool errorTrack>
    inline Scalar<1, errorTrack> operator++( //NOLINT confusing-warning
            Scalar<1, errorTrack>& s, int) {
        Scalar<1, errorTrack> temp(s);
        s += 1.0;
        return temp;
    }

    template<bool errorTrack>
    inline Scalar<1, errorTrack> operator--( //NOLINT confusing-warning
            Scalar<1, errorTrack>& s, int) {
        Scalar<1, errorTrack> temp(s);
        s -= 1.0;
        return temp;
    }
    ////////////////////////////////////////Double-WithoutError///////////////////////////////////////////
    inline Scalar<1, false>::Scalar() : d(0) {}

    inline Scalar<1, false>::Scalar(double d) : d(d) {}
    /////////////////////////////////////////Double-WithError///////////////////////////////////////////
    inline Scalar<1, true>::Scalar() : Scalar<1, false>(), a(0) {}

    inline Scalar<1, true>::Scalar(double d, double a) : Scalar<1, false>(d), a(a) {}
    //////////////////////////////////////////////Global//////////////////////////////////////////////
    template<size_t maxPrecision, bool errorTrack>
    std::ostream& operator<<(std::ostream& os, const Scalar<maxPrecision, errorTrack>& n) {
        return os << std::setprecision(10) //10 is the max precision of double.
                  << double(n)
                  << std::setprecision(6); //6 is the default precision.
    }

    template<size_t maxPrecision, bool errorTrack>
    inline Scalar<maxPrecision, errorTrack> operator+(const Scalar<maxPrecision, errorTrack>& s) {
        return Scalar<maxPrecision, errorTrack>(s);
    }

    template<size_t maxPrecision, bool errorTrack>
    inline void operator+=(Scalar<maxPrecision, errorTrack>& s1
            , const Scalar<maxPrecision, errorTrack>& s2) { s1 = s1 + s2; }

    template<size_t maxPrecision, bool errorTrack>
    inline void operator-=(Scalar<maxPrecision, errorTrack>& s1
            , const Scalar<maxPrecision, errorTrack>& s2) { s1 = s1 - s2; }

    template<size_t maxPrecision, bool errorTrack>
    inline void operator*=(Scalar<maxPrecision, errorTrack>& s1
            , const Scalar<maxPrecision, errorTrack>& s2) { s1 = s1 * s2; }

    template<size_t maxPrecision, bool errorTrack>
    inline void operator/=(Scalar<maxPrecision, errorTrack>& s1
            , const Scalar<maxPrecision, errorTrack>& s2) { s1 = s1 / s2; }

    template<size_t maxPrecision, bool errorTrack>
    inline void operator^=(Scalar<maxPrecision, errorTrack>& s1
            , const Scalar<maxPrecision, errorTrack>& s2) { s1 = s1 ^ s2; }

    template<size_t maxPrecision, bool errorTrack>
    inline void operator<<=(Scalar<maxPrecision, errorTrack>& s
            , int bits) { s = s << bits; }

    template<size_t maxPrecision, bool errorTrack>
    inline void operator>>=(Scalar<maxPrecision, errorTrack>& s
            , int bits) { s = s >> bits; }
    
    template<size_t maxPrecision, bool errorTrack>
    inline void swap(
            Scalar<maxPrecision, errorTrack>& s1
            , Scalar<maxPrecision, errorTrack>& s2
    ) noexcept {
        s1.swap(s2);
    }
}

#endif