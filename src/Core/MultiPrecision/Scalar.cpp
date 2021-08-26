/*
 * Copyright 2019-2021 WeiBo He.
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
#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Logger/LoggerRuntime.h"
#include "Physica/Core/MultiPrecision/BasicImpl/Convert.h"

namespace Physica::Core {
    namespace Internal {
        AbstractScalar<MultiPrecision>::AbstractScalar() noexcept
                : byte(nullptr), length(0), power(0) {}

        AbstractScalar<MultiPrecision>::AbstractScalar(int length_, int power_)
                : byte(reinterpret_cast<MPUnit*>(malloc(std::abs(length_) * sizeof(MPUnit))))
                , length(length_), power(power_) {
            /*
            * Length of scalar must not equal to INT_MIN or -length will make no sense.
            */
            assert(length != INT_MIN);
        }

        AbstractScalar<MultiPrecision>::AbstractScalar(const AbstractScalar<MultiPrecision>& s)
                : byte(reinterpret_cast<MPUnit*>(malloc(s.getSize() * sizeof(MPUnit))))
                , length(s.length), power(s.power) {
            memcpy(byte, s.byte, getSize() * sizeof(MPUnit));
        }

        AbstractScalar<MultiPrecision>::AbstractScalar(AbstractScalar<MultiPrecision>&& s) noexcept
                : byte(s.byte), length(s.length), power(s.power) {
            s.byte = nullptr;
        }

        AbstractScalar<MultiPrecision>::AbstractScalar(int i) : AbstractScalar(static_cast<SignedMPUnit>(i)) {}

        AbstractScalar<MultiPrecision>::AbstractScalar(SignedMPUnit unit)
                : byte(reinterpret_cast<MPUnit*>(malloc(sizeof(MPUnit))))
                , length(unit > 0 ? 1 : -1), power(0) {
            byte[0] = unit > 0 ? unit : -unit;
        }

        AbstractScalar<MultiPrecision>::AbstractScalar(double d) {
            if(d == 0) {
                byte = reinterpret_cast<MPUnit*>(malloc(sizeof(MPUnit)));
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
                byte = reinterpret_cast<MPUnit*>(malloc(length * sizeof(MPUnit)));
                //Hidden bit
                byte[1] = 1;
                byte[1] <<= remainder;
                if(remainder <= 20) {
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
                byte = reinterpret_cast<MPUnit*>(malloc(sizeof(MPUnit)));
                //Hidden bit
                byte[0] = 1;
                byte[0] <<= 20U;
                byte[0] += static_cast<MPUnit>(extract.high);
                byte[0] <<= 32U;
                byte[0] += static_cast<MPUnit>(extract.low);
                byte[0] <<= remainder - 52;
            }
        #endif
        #ifdef PHYSICA_32BIT
            if(remainder < 20) {
                length = 3;
                byte = reinterpret_cast<MPUnit*>(malloc(length * sizeof(MPUnit)));
                //Hidden bit
                byte[2] = 1;
                byte[2] <<= remainder;
                byte[2] += static_cast<MPUnit>(extract.high) >> (20 - remainder);
                byte[1] = static_cast<MPUnit>(extract.high) << (32 - (20 - remainder));
                byte[1] +=  static_cast<MPUnit>(extract.low) >> (20 - remainder);
                byte[0] = static_cast<MPUnit>(extract.low) << remainder;
            }
            else {
                length = 2;
                byte = reinterpret_cast<MPUnit*>(malloc(length * sizeof(MPUnit)));
                //Hidden bit
                byte[1] = 1;
                byte[1] <<= remainder;
                byte[1] += static_cast<MPUnit>(extract.high) << (remainder - 20);
                byte[1] += static_cast<MPUnit>(extract.low) >> (32 - (remainder - 20));
                byte[0] = static_cast<MPUnit>(extract.low) << (remainder - 20);
            }
        #endif
            if(extract.sign)
                length = -length;
        }

        AbstractScalar<MultiPrecision>::AbstractScalar(const Integer& i)
                : byte(reinterpret_cast<MPUnit*>(malloc(i.getSize() * sizeof(MPUnit))))
                , length(i.getLength())
                , power(i.getSize() - 1) {
            memcpy(byte, i.getByte(), getSize() * sizeof(MPUnit));
        }

        AbstractScalar<MultiPrecision>::AbstractScalar(const Rational& r) {
            Scalar<MultiPrecision, false> numerator(r.getNumerator());
            Scalar<MultiPrecision, false> denominator(r.getDenominator());
            Scalar<MultiPrecision, false> result = numerator / denominator;
            byte = result.byte;
            result.byte = nullptr;
            length = result.length;
            power = result.power;
        }
        /*!
        * Not accurate.
        */
        AbstractScalar<MultiPrecision>::AbstractScalar(const char* s) : AbstractScalar(strtod(s, nullptr)) {}
        /*!
        * Not accurate.
        */
        AbstractScalar<MultiPrecision>::AbstractScalar(const wchar_t* s) {
            size_t size = wcslen(s);
            char str[size + 1];
            str[size] = '\0';
            for(size_t i = 0; i < size; ++i)
                str[i] = (char)s[i];
            AbstractScalar<MultiPrecision> temp(str);
            byte = temp.byte;
            temp.byte = nullptr;
            length = temp.length;
            power = temp.power;
        }

        AbstractScalar<MultiPrecision>::~AbstractScalar() {
            free(byte);
        }

        AbstractScalar<MultiPrecision>& AbstractScalar<MultiPrecision>::operator= (
                const AbstractScalar<MultiPrecision>& s) {
            if(this == &s)
                return *this;
            length = s.length;
            int size = getSize();
            this->~AbstractScalar();
            byte = reinterpret_cast<MPUnit*>(malloc(size * sizeof(MPUnit)));
            memcpy(byte, s.byte, size * sizeof(MPUnit));
            power = s.power;
            return *this;
        }

        AbstractScalar<MultiPrecision>& AbstractScalar<MultiPrecision>::operator=(
                AbstractScalar<MultiPrecision>&& s) noexcept {
            this->~AbstractScalar();
            byte = s.byte;
            s.byte = nullptr;
            length = s.length;
            power = s.power;
            return *this;
        }

        AbstractScalar<MultiPrecision>::operator double() const {
            if(isZero())
                return 0.0;
            return Internal::convertDoubleImpl(length, power, byte);
        }

        AbstractScalar<MultiPrecision> AbstractScalar<MultiPrecision>::operator-() const {
            AbstractScalar result(-length, power);
            memcpy(result.byte, byte, getSize() * sizeof(MPUnit));
            return result;
        }

        void AbstractScalar<MultiPrecision>::swap(AbstractScalar<MultiPrecision>& s) noexcept {
            std::swap(byte, s.byte);
            std::swap(length, s.length);
            std::swap(power, s.power);
        }

        bool AbstractScalar<Float>::isInteger() const {
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

        bool AbstractScalar<Double>::isInteger() const {
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
    }
    //////////////////////////////////MultiPrecision-WithoutError///////////////////////////////////
    Scalar<MultiPrecision, false>::Scalar(const Scalar<MultiPrecision, true>& s) : Base(s) {}

    Scalar<MultiPrecision, false>::Scalar(Scalar<MultiPrecision, true>&& s) : Base(std::move(s)) {}

    Scalar<MultiPrecision, false>& Scalar<MultiPrecision, false>::operator=(const Scalar<MultiPrecision, true>& s) {
        Base::operator=(s);
        return *this;
    }

    Scalar<MultiPrecision, false>& Scalar<MultiPrecision, false>::operator=(Scalar<MultiPrecision, true>&& s) noexcept {
        Base::operator=(std::move(s));
        return *this;
    }

    Scalar<MultiPrecision, false> Scalar<MultiPrecision, false>::operator+(
            const Scalar<MultiPrecision, false>& s) const {
        auto result = addNoError(*this, s);
        cutLength(result);
        return result;
    }

    Scalar<MultiPrecision, false> Scalar<MultiPrecision, false>::operator-(
            const Scalar<MultiPrecision, false>& s) const {
        auto result = subNoError(*this, s);
        cutLength(result);
        return result;
    }

    Scalar<MultiPrecision, false> Scalar<MultiPrecision, false>::operator*(
            const Scalar<MultiPrecision, false>& s) const {
        auto result = mulNoError(*this, s);
        cutLength(result);
        return result;
    }

    Scalar<MultiPrecision, false> Scalar<MultiPrecision, false>::operator/(
            const Scalar<MultiPrecision, false>& s) const {
        return divNoError(*this, s);
    }

    Scalar<MultiPrecision, true> Scalar<MultiPrecision, false>::operator*(
            const Scalar<MultiPrecision, true>& s) const {
        auto result = mulWithError(*this, s);
        cutLength(result);
        if(s.getA() != 0)
            result.applyError(*this * s.getAccuracy());
        return result;
    }

    Scalar<MultiPrecision, true> Scalar<MultiPrecision, false>::operator/(
            const Scalar<MultiPrecision, true>& s) const {
        auto result = divWithError(*this, s);
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
        const int quotient = bits / MPUnitWidth; //NOLINT: quotient < INT_MAX
        const unsigned int remainder = bits - quotient * MPUnitWidth;
        //If remainder = 0, we must return directly because shifting a MPUnit for MPUnitWidth bits is a undefined behavior.
        if(remainder == 0) {
            Scalar copy(*this);
            copy.power += quotient;
            return copy;
        }

        const bool carry = countLeadingZeros(byte[size - 1]) < remainder;
        Scalar result(length >= 0 ? (size + carry) : -(size + carry), power + quotient + carry);
        result.byte[0] = 0;
        const int size_1 = size - 1;
        for(int i = 0; i < size_1; ++i) {
            result.byte[i] |= byte[i] << remainder;
            result.byte[i + 1] = byte[i] >> (MPUnitWidth - remainder);
        }
        result.byte[size_1] |= byte[size_1] << remainder;
        if(carry)
            result.byte[size] = byte[size_1] >> (MPUnitWidth - remainder);
        return result;
    }

    Scalar<MultiPrecision, false> Scalar<MultiPrecision, false>::operator>>(int bits) const {
        if(bits == 0)
            return Scalar(*this);
        if(bits < 0)
            return *this << -bits;
        const int size = getSize();
        const int quotient = bits / MPUnitWidth; //NOLINT: quotient < INT_MAX
        const unsigned int remainder = bits - quotient * MPUnitWidth;
        //If remainder = 0, we must return directly because shifting a MPUnit for MPUnitWidth bits is a undefined behavior.
        if(remainder == 0) {
            Scalar copy(*this);
            copy.power -= quotient;
            return copy;
        }

        const bool carry = (countLeadingZeros(byte[size - 1]) + remainder) < MPUnitWidth;
        Scalar result(length >= 0 ? (size + carry) : -(size + carry), power - quotient + carry - 1);
        if(carry)
            result.byte[size] = byte[size - 1] >> remainder;

        for(int i = size - 1; i > 0; --i) {
            result.byte[i] = byte[i] << (MPUnitWidth - remainder);
            result.byte[i] |= byte[i - 1] >> remainder;
        }
        result.byte[0] = byte[0] << (MPUnitWidth - remainder);
        return result;
    }

    Scalar<MultiPrecision, false> Scalar<MultiPrecision, false>::operator-() const {
        return static_cast<Scalar<MultiPrecision, false>&&>(Base::operator-());
    }
    /*!
     * return true if the abstract value of s1 is larger or equal than the abstract value of s2.
     * return false if the abstract value of s1 is smaller to the abstract value of s2.
     *
     * Optimize:
     * Is subtract faster than comparing the elements?
     */
    bool absCompare(const Internal::AbstractScalar<MultiPrecision>& s1, const Internal::AbstractScalar<MultiPrecision>& s2) {
        if(s1.isZero() || s2.isZero())
            return true;
        if(s1.getPower() > s2.getPower())
            return true;
        if(s1.getPower() < s2.getPower())
            return false;
        const Internal::AbstractScalar<MultiPrecision>* longer, *shorter;
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
    bool operator>(const Internal::AbstractScalar<MultiPrecision>& s1, const Internal::AbstractScalar<MultiPrecision>& s2) {
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
            auto scalar1 = static_cast<const Scalar<MultiPrecision, false>&>(s1);
            auto scalar2 = static_cast<const Scalar<MultiPrecision, false>&>(s2);
            //The only method left.
            //Optimize: We have confirmed that s1, s2 have the same sign and power, possible make use them to get better performance.
            Scalar<MultiPrecision, false> subtract = scalar1 - scalar2;
            result = subtract.isPositive();
        }
        return result;
    }
    /*!
     * Optimize:
     * Is subtract faster than comparing the elements?
     */
    bool operator<(const Internal::AbstractScalar<MultiPrecision>& s1, const Internal::AbstractScalar<MultiPrecision>& s2) {
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
            auto scalar1 = static_cast<const Scalar<MultiPrecision, false>&>(s1);
            auto scalar2 = static_cast<const Scalar<MultiPrecision, false>&>(s2);
            //The only method left.
            //Optimize: We have confirmed that s1, s2 have the same sign and power, possible make use them to get better performance.
            Scalar<MultiPrecision, false> subtract = scalar1 - scalar2;
            result = subtract.isNegative();
        }
        return result;
    }

    bool operator==(const Internal::AbstractScalar<MultiPrecision>& s1, const Internal::AbstractScalar<MultiPrecision>& s2) {
        auto scalar1 = static_cast<const Scalar<MultiPrecision, false>&>(s1);
        auto scalar2 = static_cast<const Scalar<MultiPrecision, false>&>(s2);
        return s1.getPower() == s2.getPower()
               //Here length may not equal n.length because we define numbers like 1.0 and 1.00 are equal.
               && ((s1.getLength() ^ s2.getLength()) >= 0) //NOLINT
               //Optimize: We have confirmed that s1, s2 have the same sign and power, possible make use them to get better performance.
               && Scalar<MultiPrecision, false>(scalar1 - scalar2).isZero();
    }
    ///////////////////////////////////MultiPrecision-WithError/////////////////////////////////////
    Scalar<MultiPrecision, true>::Scalar() noexcept : AbstractScalar<MultiPrecision>(), a(0) {}

    Scalar<MultiPrecision, true>::Scalar(int length_, int power_, MPUnit a_) noexcept
            : AbstractScalar<MultiPrecision>(length_, power_), a(a_) {}

    Scalar<MultiPrecision, true>::Scalar(const Scalar<MultiPrecision, true>& s) : AbstractScalar<MultiPrecision>(s) {
        a = s.a;
    }

    Scalar<MultiPrecision, true>::Scalar(Scalar<MultiPrecision, true>&& s) noexcept
            : AbstractScalar<MultiPrecision>(std::move(s)), a(s.a) {}

    Scalar<MultiPrecision, true>::Scalar(int i, MPUnit a_)
            : AbstractScalar<MultiPrecision>(i), a(a_) {}

    Scalar<MultiPrecision, true>::Scalar(SignedMPUnit unit, MPUnit a_)
            : AbstractScalar<MultiPrecision>(unit), a(a_) {}

    Scalar<MultiPrecision, true>::Scalar(double d, MPUnit a_)
            : AbstractScalar<MultiPrecision>(d), a(a_) {}

    Scalar<MultiPrecision, true>::Scalar(const char* s, MPUnit a_)
            : AbstractScalar<MultiPrecision>(s), a(a_) {}

    Scalar<MultiPrecision, true>::Scalar(const wchar_t* s, MPUnit a_)
            : AbstractScalar<MultiPrecision>(s), a(a_) {}

    Scalar<MultiPrecision, true>::Scalar(const Scalar<MultiPrecision, false>& s, MPUnit a_) : Base(s), a(a_) {}

    Scalar<MultiPrecision, true>::Scalar(Scalar<MultiPrecision, false>&& s, MPUnit a_) : Base(std::move(s)), a(a_) {}

    Scalar<MultiPrecision, true>& Scalar<MultiPrecision, true>::operator=(const Scalar<MultiPrecision, false>& s) {
        Base::operator=(s);
        a = 0;
        return *this;
    }

    Scalar<MultiPrecision, true>& Scalar<MultiPrecision, true>::operator=(Scalar<MultiPrecision, false>&& s) noexcept {
        Base::operator=(std::move(s));
        a = 0;
        return *this;
    }

    Scalar<MultiPrecision, true> Scalar<MultiPrecision, true>::operator*(
            const Scalar<MultiPrecision, false>& s) const {
        Q_UNUSED(MultiPrecision)
        auto result = mulWithError(*this, s);
        cutLength(result);
        if(getA() != 0)
            result.applyError(s * getAccuracy());
        return result;
    }

    Scalar<MultiPrecision, true> Scalar<MultiPrecision, true>::operator/(
            const Scalar<MultiPrecision, false>& s) const {
        Q_UNUSED(MultiPrecision)
        auto result = divWithError(*this, s);
        if(getA() != 0)
            result.applyError(getAccuracy() / s);
        return result;
    }

    Scalar<MultiPrecision, true> Scalar<MultiPrecision, true>::operator+(
            const Scalar<MultiPrecision, true>& s) const {
        auto result = addWithError(*this, s);
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
        auto result = subWithError(*this, s);
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
        auto result = mulWithError(*this, s);
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
        auto result = divWithError(*this, s);
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
        int size = getSize();
        const int quotient = bits / MPUnitWidth; //NOLINT: quotient < INT_MAX
        const unsigned int remainder = bits - quotient * MPUnitWidth;
        //If remainder = 0, we must return directly because shifting a MPUnit for MPUnitWidth bits is a undefined behavior.
        if(remainder == 0) {
            Scalar copy(*this);
            copy.power += quotient;
            return copy;
        }

        const int size_1 = size - 1;
        const bool carry = countLeadingZeros(byte[size_1]) < remainder;
        const bool a_carry = countLeadingZeros(a) < remainder;
        Scalar result(length >= 0 ? (size + carry - a_carry) : -(size + carry - a_carry), power + quotient + carry,
                      a_carry ? a >> (MPUnitWidth - remainder) : a << remainder); //Add 1 to a if a_carry is true to get more accurate estimation.
        const auto byte_start = byte + a_carry;
        size -= a_carry;
        result.byte[0] = 0;
        for(int i = 0; i < size_1; ++i) {
            result.byte[i] |= byte_start[i] << remainder;
            result.byte[i + 1] = byte_start[i] >> (MPUnitWidth - remainder);
        }
        result.byte[size_1] |= byte_start[size_1] << remainder;
        if(carry)
            result.byte[size] = byte_start[size_1] >> (MPUnitWidth - remainder);
        return result;
    }

    Scalar<MultiPrecision, true> Scalar<MultiPrecision, true>::operator>>(int bits) const {
        if(bits == 0)
            return Scalar(*this);
        if(bits < 0)
            return *this << -bits;
        const int size = getSize();
        const int quotient = bits / MPUnitWidth; //NOLINT: quotient < INT_MAX
        const unsigned int remainder = bits - quotient * MPUnitWidth;
        //If remainder = 0, we must return directly because shifting a MPUnit for MPUnitWidth bits is a undefined behavior.
        if(remainder == 0) {
            Scalar copy(*this);
            copy.power -= quotient;
            return copy;
        }

        const int size_1 = size - 1;
        const bool carry = (countLeadingZeros(byte[size_1]) + remainder) < MPUnitWidth;
        const MPUnit accuracy = a >> remainder;
        Scalar result(length >= 0 ? (size + carry) : -(size + carry)
                , power - quotient + carry - 1, accuracy > 0 ? accuracy : 0);
        if(carry)
            result.byte[size] = byte[size_1] >> remainder;

        for(int i = size_1; i > 0; --i) {
            result.byte[i] = byte[i] << (MPUnitWidth - remainder);
            result.byte[i] |= byte[i - 1] >> remainder;
        }
        result.byte[0] = byte[0] << (MPUnitWidth - remainder);
        return result;
    }

    Scalar<MultiPrecision, true> Scalar<MultiPrecision, true>::operator-() const {
        return Scalar<MultiPrecision, true>(static_cast<Scalar<MultiPrecision, false>&&>(Base::operator-()), a);
    }
    /*!
     * Add error to this and adjust this->length as well as this->byte.
     */
    Scalar<MultiPrecision, true>& Scalar<MultiPrecision, true>::applyError(
            const Scalar<MultiPrecision, false>& error) {
        assert(!error.isNegative());
        if(power == error.getPower() || (power > error.getPower() && power - error.getPower() > 0)) {
            int size = getSize();
            const int copy = size;
            const int temp = power - error.getPower() - size + 1;//Equals to (power - size + 1 - error.getPower()). Exchanged the order to avoid overflow.
            MPUnit copy_a = a;
            if(temp <= 0) {
                if(temp < 0) {
                    Scalar<MultiPrecision, false> error_1 = getAccuracy() + error;
                    size += temp;
                    //Use (a += error_1.byte[error_1.getSize() - 1] + 1) for more conservative error estimate.
                    a += error_1[error_1.getSize() - 1];
                }
                else
                    //Use (a += error.byte[error.getSize() - 1] + 1) for more conservative error estimate.
                    a += error[error.getSize() - 1];
            }

            if(a < copy_a) {
                //Use (a = 2) for more conservative error estimate.
                a = 1;
                --size;
            }

            if(size < 1) {
                power += 1 - size;
                goto errorOverflow;
            }

            if(size < copy) {
                auto new_byte = reinterpret_cast<MPUnit*>(malloc(size * sizeof(MPUnit)));
                memcpy(new_byte, byte + copy - size, size * sizeof(MPUnit));
                this->~Scalar();
                byte = new_byte;
            }
            length = length < 0 ? -size : size;
        }
        else if(power < error.getPower()) {
            power = error.getPower();
            goto errorOverflow;
        }
        return *this;
    //Handle error overflow, power must be set before goto here.
    errorOverflow:
        Warning(0, "Accumulated too many errors.");
        length = length < 0 ? -1 : 1;
        return *this;
    }

    void Scalar<MultiPrecision, true>::swap(
            Scalar<MultiPrecision, true>& s) noexcept {
        Base::swap(s);
        std::swap(a, s.a);
    }
    //!Return accuracy in class Scalar.
    Scalar<MultiPrecision, false> Scalar<MultiPrecision, true>::getAccuracy() const {
        Scalar<MultiPrecision, false> result(1, power - Scalar<MultiPrecision, true>::getSize() + 1);
        result.byte[0] = a;
        return result;
    }
    ///////////////////////////////////////////Float-Double////////////////////////////////////////////////
    void Scalar<Float, true>::swap(Scalar<Float, true>& s) noexcept {
        Base::swap(s);
        std::swap(a, s.a);
    }

    void Scalar<Double, true>::swap(Scalar<Double, true>& s) noexcept {
        Base::swap(s);
        std::swap(a, s.a);
    }
    ///////////////////////////////////////////Global////////////////////////////////////////////////
    template<>
    std::ostream& operator<<(std::ostream& os, const Scalar<MultiPrecision, false>& s) {
        const auto& basicConst = BasicConst::getInstance();
        const int power = s.getPower();
        int exp = int(power * basicConst.ln_2_10);
        double coefficient = std::exp(power * basicConst.ln_2 - exp * basicConst.ln_10) * s[s.getSize() - 1];
        while(coefficient > 10) {
            ++exp;
            coefficient /= 10;
        }
        os << coefficient << " * 10 ^ " << exp;
        return os;
    }
}