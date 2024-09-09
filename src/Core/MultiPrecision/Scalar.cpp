/*
 * Copyright 2019-2024 Weibo He.
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
    //////////////////////////////////MultiPrecision///////////////////////////////////
    Scalar<MultiPrecision>::Scalar() : byte(nullptr), length(0), power(0) {}

    Scalar<MultiPrecision>::Scalar(int length_, int power_)
            : byte(reinterpret_cast<MPUnit*>(malloc(std::abs(length_) * sizeof(MPUnit))))
            , length(length_), power(power_) {
        /**
         * Length of scalar must not equal to INT_MIN or -length will make no sense.
         */
        assert(length != INT_MIN);
    }

    Scalar<MultiPrecision>::Scalar(int i) : Scalar(static_cast<SignedMPUnit>(i)) {}

    Scalar<MultiPrecision>::Scalar(SignedMPUnit unit)
            : byte(reinterpret_cast<MPUnit*>(malloc(sizeof(MPUnit))))
            , length(unit > 0 ? 1 : -1), power(0) {
            byte[0] = unit > 0 ? unit : -unit;
    }

    Scalar<MultiPrecision>::Scalar(double d) {
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
        if constexpr (PhysicaWordSize == 64) {
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
        }
        else {
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
        }

        if(extract.sign)
            length = -length;
    }

    Scalar<MultiPrecision>::Scalar(const Integer& i)
            : byte(reinterpret_cast<MPUnit*>(malloc(i.getSize() * sizeof(MPUnit))))
            , length(i.getLength())
            , power(i.getSize() - 1) {
        memcpy(byte, i.getByte(), getSize() * sizeof(MPUnit));
    }

    Scalar<MultiPrecision>::Scalar(const Rational& r) {
        Scalar<MultiPrecision> numerator(r.getNumerator());
        Scalar<MultiPrecision> denominator(r.getDenominator());
        Scalar<MultiPrecision> result = numerator / denominator;
        byte = result.byte;
        result.byte = nullptr;
        length = result.length;
        power = result.power;
    }
    /**
     * Not accurate.
     */
    Scalar<MultiPrecision>::Scalar(const char* s) : Scalar(strtod(s, nullptr)) {}
    /**
    * Not accurate.
    */
    Scalar<MultiPrecision>::Scalar(const wchar_t* s) {
        size_t size = wcslen(s);
        char str[size + 1];
        str[size] = '\0';
        for(size_t i = 0; i < size; ++i)
            str[i] = (char)s[i];
        Scalar<MultiPrecision> temp(str);
        byte = temp.byte;
        temp.byte = nullptr;
        length = temp.length;
        power = temp.power;
    }

    MPUnit Scalar<MultiPrecision>::operator[](unsigned int index) const {
        assert(index < static_cast<unsigned int>(getSize()));
        return byte[index];
    }

    Scalar<MultiPrecision>::operator double() const {
        if(isZero())
            return 0.0;
        return Internal::convertDoubleImpl(length, power, byte);
    }

    Scalar<MultiPrecision>::Scalar(const Scalar<MultiPrecision>& s)
            : byte(reinterpret_cast<MPUnit*>(malloc(s.getSize() * sizeof(MPUnit))))
            , length(s.length), power(s.power) {
        memcpy(byte, s.byte, getSize() * sizeof(MPUnit));
    }

    Scalar<MultiPrecision>::Scalar(Scalar<MultiPrecision>&& s) noexcept
            : byte(s.byte), length(s.length), power(s.power) {
        s.byte = nullptr;
    }

    Scalar<MultiPrecision>::~Scalar() {
        free(byte);
    }

    Scalar<MultiPrecision> Scalar<MultiPrecision>::operator+(
            const Scalar<MultiPrecision>& s) const {
        auto result = add(*this, s);
        cutLength(result);
        return result;
    }

    Scalar<MultiPrecision> Scalar<MultiPrecision>::operator-(
            const Scalar<MultiPrecision>& s) const {
        auto result = sub(*this, s);
        cutLength(result);
        return result;
    }

    Scalar<MultiPrecision> Scalar<MultiPrecision>::operator*(
            const Scalar<MultiPrecision>& s) const {
        auto result = mul(*this, s);
        cutLength(result);
        return result;
    }

    Scalar<MultiPrecision> Scalar<MultiPrecision>::operator/(
            const Scalar<MultiPrecision>& s) const {
        return div(*this, s);
    }

    Scalar<MultiPrecision> Scalar<MultiPrecision>::operator<<(int bits) const {
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

    Scalar<MultiPrecision> Scalar<MultiPrecision>::operator>>(int bits) const {
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

    Scalar<MultiPrecision> Scalar<MultiPrecision>::operator-() const {
        Scalar result(-length, power);
        memcpy(result.byte, byte, getSize() * sizeof(MPUnit));
        return result;
    }

    void Scalar<MultiPrecision>::swap(Scalar<MultiPrecision>& __restrict s) noexcept {
        assert(this != &s && "[Error]: Self swap is likely a bug");
        std::swap(byte, s.byte);
        std::swap(length, s.length);
        std::swap(power, s.power);
    }
    /*!
     * return true if the abstract value of s1 is larger or equal than the abstract value of s2.
     * return false if the abstract value of s1 is smaller to the abstract value of s2.
     *
     * Optimize:
     * Is subtract faster than comparing the elements?
     */
    bool absCompare(const Scalar<MultiPrecision>& s1, const Scalar<MultiPrecision>& s2) {
        if(s1.isZero() || s2.isZero())
            return true;
        if(s1.getPower() > s2.getPower())
            return true;
        if(s1.getPower() < s2.getPower())
            return false;
        const Scalar<MultiPrecision>* longer, *shorter;
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
    bool operator>(const Scalar<MultiPrecision>& s1, const Scalar<MultiPrecision>& s2) {
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
            auto scalar1 = static_cast<const Scalar<MultiPrecision>&>(s1);
            auto scalar2 = static_cast<const Scalar<MultiPrecision>&>(s2);
            //The only method left.
            //Optimize: We have confirmed that s1, s2 have the same sign and power, possible make use them to get better performance.
            Scalar<MultiPrecision> subtract = scalar1 - scalar2;
            result = subtract.isPositive();
        }
        return result;
    }
    /*!
     * Optimize:
     * Is subtract faster than comparing the elements?
     */
    bool operator<(const Scalar<MultiPrecision>& s1, const Scalar<MultiPrecision>& s2) {
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
            auto scalar1 = static_cast<const Scalar<MultiPrecision>&>(s1);
            auto scalar2 = static_cast<const Scalar<MultiPrecision>&>(s2);
            //The only method left.
            //Optimize: We have confirmed that s1, s2 have the same sign and power, possible make use them to get better performance.
            Scalar<MultiPrecision> subtract = scalar1 - scalar2;
            result = subtract.isNegative();
        }
        return result;
    }

    bool operator==(const Scalar<MultiPrecision>& s1, const Scalar<MultiPrecision>& s2) {
        auto scalar1 = static_cast<const Scalar<MultiPrecision>&>(s1);
        auto scalar2 = static_cast<const Scalar<MultiPrecision>&>(s2);
        return s1.getPower() == s2.getPower()
               //Here length may not equal n.length because we define numbers like 1.0 and 1.00 are equal.
               && ((s1.getLength() ^ s2.getLength()) >= 0) //NOLINT
               //Optimize: We have confirmed that s1, s2 have the same sign and power, possible make use them to get better performance.
               && Scalar<MultiPrecision>(scalar1 - scalar2).isZero();
    }
    ///////////////////////////////////////////Float////////////////////////////////////////////////
    std::istream& operator>>(std::istream& is, Scalar<Float>& scalar) {
        is >> scalar.f;
        return is;
    }

    bool Scalar<Float>::isInteger() const {
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
        /**
         * exp + zeros - 127 >= 23
         * , 127 is the exp bias of float numbers, 23 is the number of bits of significand of float numbers.
         * We move 127 to the right hand side to avoid underflow.
         */
        return exp + zeros >= 150;
    }
    ///////////////////////////////////////////Double////////////////////////////////////////////////
    std::istream& operator>>(std::istream& is, Scalar<Double>& scalar) {
        is >> scalar.d;
        return is;
    }

    bool Scalar<Double>::isInteger() const {
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
        /**
         * exp + zeros - 1023 >= 52
         * , 1023 is the exp bias of float numbers, 52 is the number of bits of significand of float numbers.
         * We move 1023 to the right hand side to avoid underflow.
         */
        return exp + zeros >= 1075;
    }
    ///////////////////////////////////////////Global////////////////////////////////////////////////
    std::ostream& operator<<(std::ostream& os, const Scalar<MultiPrecision>& s) {
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