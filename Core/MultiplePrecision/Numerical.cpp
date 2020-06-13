/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <cmath>
#include <iostream>
#include <iomanip>
#include "Operation/Pow.h"

namespace Physica::Core {
    ////////////////////////////////Numerical////////////////////////////////
    Numerical::Numerical() noexcept : byte(nullptr), length(0), power(0), a(0) {}

    Numerical::Numerical(int length, int power, NumericalUnit a) noexcept
            : byte(reinterpret_cast<NumericalUnit*>(malloc(abs(length) * sizeof(NumericalUnit))))
            , length(length), power(power), a(a) {}
    /*
     * byte must be allocated by malloc().
     * It is declared to be protected to avoid incorrect use.
     */
    Numerical::Numerical(NumericalUnit* byte, int length, int power, NumericalUnit a)
            : byte(byte), length(length), power(power), a(a) {}

    Numerical::Numerical(const Numerical& n) noexcept
            : byte(reinterpret_cast<NumericalUnit*>(malloc(n.getSize() * sizeof(NumericalUnit))))
            , length(n.length), power(n.power), a(n.a) {
        memcpy(byte, n.byte, getSize() * sizeof(NumericalUnit));
    }

    Numerical::Numerical(Numerical&& n) noexcept : byte(n.byte), length(n.length), power(n.power), a(n.a) {
        n.byte = nullptr;
    }

    Numerical::Numerical(SignedNumericalUnit unit, NumericalUnit a) noexcept
            : byte(reinterpret_cast<NumericalUnit*>(malloc(sizeof(NumericalUnit)))), power(0), a(a) {
        if(unit < 0) {
            byte[0] = -unit;
            length = -1;
        }
        else {
            byte[0] = unit;
            length = 1;
        }
    }

    Numerical::Numerical(double d, NumericalUnit a) : a(a) {
        if(d == 0) {
            byte = reinterpret_cast<NumericalUnit*>(malloc(sizeof(NumericalUnit)));
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
            byte = reinterpret_cast<NumericalUnit*>(malloc(length * sizeof(NumericalUnit)));
            //Hidden bit
            byte[1] = 1;
            byte[1] <<= remainder;
            if(remainder <= 20) {
                byte[1] += static_cast<NumericalUnit>(extract.structure.high) >> (20 - remainder);
                byte[0] = static_cast<NumericalUnit>(extract.structure.high) << (44 + remainder);
                byte[0] += static_cast<NumericalUnit>(extract.structure.low) << (32 - (20 - remainder));
            }
            else {
                byte[1] += static_cast<NumericalUnit>(extract.structure.high) << (remainder - 20);
                byte[1] += static_cast<NumericalUnit>(extract.structure.low) >> (32 - (remainder - 20));
                byte[0] = static_cast<NumericalUnit>(extract.structure.low) << (32 + (remainder - 20));
            }
        }
        else {
            length = 1;
            byte = reinterpret_cast<NumericalUnit*>(malloc(sizeof(NumericalUnit)));
            //Hidden bit
            byte[0] = 1;
            byte[0] <<= 20U;
            byte[0] += static_cast<NumericalUnit>(extract.structure.high);
            byte[0] <<= 32U;
            byte[0] += static_cast<NumericalUnit>(extract.structure.low);
            byte[0] <<= remainder - 52;
        }
    #endif
    #ifdef PHYSICA_32BIT
        if(remainder < 20) {
        length = 3;
        byte = reinterpret_cast<NumericalUnit*>(malloc(length * sizeof(NumericalUnit)));
        //Hidden bit
        byte[2] = 1;
        byte[2] <<= remainder;
        byte[2] += static_cast<NumericalUnit>(extract.structure.high) >> (20 - remainder);
        byte[1] = static_cast<NumericalUnit>(extract.structure.high) << (32 - (20 - remainder));
        byte[1] +=  static_cast<NumericalUnit>(extract.structure.low) >> (20 - remainder);
        byte[0] = static_cast<NumericalUnit>(extract.structure.low) << remainder;
    }
    else {
        length = 2;
        byte = reinterpret_cast<NumericalUnit*>(malloc(length * sizeof(NumericalUnit)));
        //Hidden bit
        byte[1] = 1;
        byte[1] <<= remainder;
        byte[1] += static_cast<NumericalUnit>(extract.structure.high) << (remainder - 20);
        byte[1] += static_cast<NumericalUnit>(extract.structure.low) >> (32 - (remainder - 20));
        byte[0] = static_cast<NumericalUnit>(extract.structure.low) << (remainder - 20);
    }
    #endif
        if(extract.structure.sign)
            length = -length;
    }
    /*
     * Not accurate.
     */
    Numerical::Numerical(const char* s, NumericalUnit a) : Numerical(strtod(s, nullptr), a) {}
    /*
     * Not accurate.
     */
    Numerical::Numerical(const wchar_t* s, NumericalUnit a) : a(a) {
        size_t size = wcslen(s);
        char str[size + 1];
        str[size] = '\0';
        for(int i = 0; i < size; ++i)
            str[i] = (char)s[i];
        Numerical temp(str);
        byte = temp.byte;
        temp.byte = nullptr;
        length = temp.length;
        power = temp.power;
    }
    /*
     * Not accurate.
     */
    Numerical::Numerical(const std::string& s, NumericalUnit a) : Numerical(s.c_str(), a) {}
    /*
     * Not accurate.
     */
    Numerical::Numerical(const std::wstring& s, NumericalUnit a) : Numerical(s.c_str(), a) {}

    Numerical::~Numerical() { free(byte); }

    Numerical::operator double() const {
        if(isZero())
            return 0.0;
        double_extract extract{0};
        extract.structure.sign = length < 0;

        const auto zeroCount = countLeadingZeros(byte[getSize() - 1]);
        auto exp = power * __WORDSIZE + NumericalUnitWidth - zeroCount - 1 + 1023;
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

    std::ostream& operator<<(std::ostream& os, const Numerical& n) {
        //10 is the max precision of double.
        return os << std::setprecision(10) << std::to_string(double(n)) << "\tLength = "
                  << n.getSize() << "\tPower = " << n.power << "\tAccuracy = " << (int)n.a << std::setprecision(6); //6 is the default precision.
    }

    Numerical Numerical::operator+(const Numerical& n) const {
        Numerical result = add(*this, n);
        cutLength(result);
        if(getA() != 0 || n.getA() != 0) {
            if(getA() == 0)
                result.applyError(n.getAccuracy());
            else if(n.getA() == 0)
                result.applyError(getAccuracy());
            else
                result.applyError(add(getAccuracy(), n.getAccuracy()));
        }
        return result;
    }

    Numerical Numerical::operator-(const Numerical& n) const {
        Numerical result = sub(*this, n);
        cutLength(result);
        if(getA() != 0 || n.getA() != 0) {
            if(getA() == 0)
                result.applyError(n.getAccuracy());
            else if(n.getA() == 0)
                result.applyError(getAccuracy());
            else
                result.applyError(add(getAccuracy(), n.getAccuracy()));
        }
        return result;
    }

    Numerical Numerical::operator*(const Numerical& n) const {
        Numerical result = mul(*this, n);
        cutLength(result);
        if(getA() != 0 || n.getA() != 0) {
            if(getA() == 0)
                result.applyError(mul(*this, n.getAccuracy()));
            else if(n.getA() == 0)
                result.applyError(mul(n, getAccuracy()));
            else {
                Numerical n1_a = getAccuracy();
                Numerical n_a = n.getAccuracy();
                Numerical temp1 = mul(n1_a, n_a);
                Numerical temp2 = add(mul(*this, n_a), mul(n, n1_a));
                result.applyError(add(temp1, temp2));
            }
        }
        return result;
    }

    Numerical Numerical::operator/(const Numerical& n) const {
        Numerical result = div(*this, n);
        if(getA() != 0 || n.getA() != 0) {
            if(getA() == 0) {
                Numerical n_a = n.getAccuracy();
                Numerical temp_1 = mul(*this, n_a);
                Numerical temp_2 = mul(n, sub(n, n_a));
                result.applyError(div(temp_1, temp_2));
            }
            else if(n.getA() == 0)
                result.applyError(div(getAccuracy(), n));
            else {
                Numerical n_a = n.getAccuracy();
                Numerical temp_1 = add(mul(*this, n_a), mul(n, getAccuracy()));
                Numerical temp_2 = mul(n, sub(n, n_a));
                result.applyError(div(temp_1, temp_2));
            }
        }
        return result;
    }
    //FIXME a is not accurate
    Numerical Numerical::operator<<(int bits) const {
        if(bits == 0)
            return Numerical(*this);
        if(bits < 0)
            return *this >> -bits;
        const int size = getSize();
        const int quotient = bits / NumericalUnitWidth; //NOLINT: quotient < INT_MAX
        const unsigned int remainder = bits - quotient * NumericalUnitWidth;
        const bool carry = countLeadingZeros(byte[size - 1]) < remainder;

        Numerical result(length >= 0 ? (size + carry) : -(size + carry), power + quotient + carry, a);
        result[0] = 0;
        for(int i = 0; i < size - 1; ++i) {
            result[i] |= byte[i] << remainder;
            result[i + 1] = byte[i] >> (NumericalUnitWidth - remainder);
        }
        result[size - 1] |= byte[size - 1] << remainder;
        if(carry)
            result[size] = byte[size - 1] >> (NumericalUnitWidth - remainder);

        return result;
    }
    //FIXME a is not accurate
    Numerical Numerical::operator>>(int bits) const {
        if(bits == 0)
            return Numerical(*this);
        if(bits < 0)
            return *this << -bits;
        const int size = getSize();
        const int quotient = bits / NumericalUnitWidth; //NOLINT: quotient < INT_MAX
        const unsigned int remainder = bits - quotient * NumericalUnitWidth;
        const bool carry = (countLeadingZeros(byte[size - 1]) + remainder) < NumericalUnitWidth;

        Numerical result(length >= 0 ? (size + carry) : -(size + carry), power - quotient + carry - 1, a);
        if(carry)
            result[size] = byte[size - 1] >> remainder;

        for(int i = size - 1; i > 0; --i) {
            result[i] = byte[i] << (NumericalUnitWidth - remainder);
            result[i] |= byte[i - 1] >> remainder;
        }
        result[0] = byte[0] << (NumericalUnitWidth - remainder);

        return result;
    }

    Numerical& Numerical::operator= (const Numerical& n) {
        if(this == &n)
            return *this;
        length = n.length;
        int size = getSize();
        this->~Numerical();
        byte = reinterpret_cast<NumericalUnit*>(malloc(size * sizeof(NumericalUnit)));
        memcpy(byte, n.byte, size * sizeof(NumericalUnit));
        power = n.power;
        a = n.a;
        return *this;
    }

    Numerical& Numerical::operator=(Numerical&& n) noexcept {
        this->~Numerical();
        byte = n.byte;
        n.byte = nullptr;
        length = n.length;
        power = n.power;
        a = n.a;
        return *this;
    }
    //Reference: MaTHmu Project Group.计算机代数系统的数学原理[M].Beijing: TsingHua University Press, 2009.45
    Numerical Numerical::operator^ (const Numerical& n) const {
        if(isInteger())
            return powNumerical(*this, n);
        return exp(ln(*this) * n);
    }

    Numerical Numerical::operator-() const {
        Numerical result(-length, power, a);
        memcpy(result.byte, byte, getSize() * sizeof(NumericalUnit));
        return result;
    }
    /*!
     * return true if the abstract value of n1 is larger or equal than the abstract value of n2.
     * return false if the abstract value of n1 is smaller to the abstract value of n2.
     *
     * Optimize:
     * Is subtract faster than comparing the elements?
     */
    bool absCompare(const Numerical& n1, const Numerical& n2) {
        if(n1.getPower() > n2.getPower())
            return true;
        if(n1.getPower() < n2.getPower())
            return false;
        const Numerical* longer, *shorter;
        int longer_length, shorter_length;
        /* Compare length */ {
            const auto n1_length = n1.getLength(), n2_length = n2.getLength();
            bool b = n1_length > n2_length;
            longer = b ? &n1 : &n2;
            shorter = b ? &n2 : &n1;
            longer_length = b ? n1_length : n2_length;
            shorter_length = b ? n2_length : n1_length;
        }
        --longer_length;
        --shorter_length;
        for(; shorter_length >= 0; --shorter_length, --longer_length) {
            if((*longer)[longer_length] > (*shorter)[shorter_length])
                return longer == &n1;
            if((*longer)[longer_length] < (*shorter)[shorter_length])
                return shorter == &n1;
        }
        return longer_length == shorter_length || longer == &n1;
    }
    /*!
     * Optimize:
     * Is subtract faster than comparing the elements?
     */
    bool operator> (const Numerical& n1, const Numerical& n2) {
        //Judge from sign.
        if(n1.isPositive()) {
            if(!n2.isPositive())
                return true;
        }
        else if(n1.isZero())
            return n2.isNegative();
        else {
            if(!n2.isNegative())
                return false;
        }
        //If we cannot get a result, judge from power
        bool result;
        if(n1.getPower() > n2.getPower())
            result = true;
        else if(n1.getPower() < n2.getPower())
            result = false;
        else {
            //The only method left.
            Numerical subtract = n1 - n2;
            result = subtract.isPositive();
        }
        return result;
    }
    /*!
     * Optimize:
     * Is subtract faster than comparing the elements?
     */
    bool operator< (const Numerical& n1, const Numerical& n2) {
        //Judge from sign.
        if(n1.isPositive()) {
            if(!n2.isPositive())
                return false;
        }
        else if(n1.isZero())
            return n2.isPositive();
        else {
            if(!n2.isNegative())
                return true;
        }
        //If we cannot get a result, judge from power
        bool result;
        if(n1.getPower() > n2.getPower())
            result = false;
        else if(n1.getPower() < n2.getPower())
            result = true;
        else {
            //The only method left.
            Numerical subtract = n1 - n2;
            result = subtract.isNegative();
        }
        return result;
    }

    bool operator>= (const Numerical& n1, const Numerical& n2) {
        return !(n1 < n2);
    }

    bool operator<= (const Numerical& n1, const Numerical& n2) {
        return !(n1 > n2);
    }

    bool operator== (const Numerical& n1, const Numerical& n2) {
        return n1.getPower() == n2.getPower()
        //Here length may not equal n.length because we define numbers like 1.0 and 1.00 are equal.
        && ((n1.getLength() ^ n2.getLength()) >= 0) //NOLINT
        && n1.getA() == n2.getA()
        && (n1 - n2).isZero();
    }

    bool operator!= (const Numerical& n1, const Numerical& n2) {
        return !(n1 == n2);
    }
    //Return accuracy in class Numerical.
    Numerical Numerical::getAccuracy() const {
        Numerical result(1, power - getSize() + 1);
        result[0] = a;
        return result;
    }
    /*
     * Add error to this and adjust this->length as well as this-> byte.
     *
     * Optimize:
     * error is always non-zero, if(!error.isZero()) is unnecessary.
     */
    Numerical& Numerical::applyError(const Numerical& error) {
        if(!error.isZero()) {
            int size = getSize();
            int copy = size;
            int temp = power - size + 1 - error.power;
            NumericalUnit copy_a = a;
            if(temp <= 0) {
                if(temp < 0) {
                    auto error_1 = add(getAccuracy(), error);
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
                auto new_byte = reinterpret_cast<NumericalUnit*>(malloc(size * sizeof(NumericalUnit)));
                memcpy(new_byte, byte + copy - size, size * sizeof(NumericalUnit));
                this->~Numerical();
                byte = new_byte;
            }
            length = length < 0 ? -size : size;
        }
        return *this;
    }

    void Numerical::swap(Numerical &n) noexcept {
        auto temp_byte = byte;
        byte = n.byte;
        n.byte = temp_byte;
        auto temp_length = length;
        length = n.length;
        n.length = temp_length;
        temp_length = power;
        power = n.power;
        n.power = temp_length;
        NumericalUnit temp_a = a;
        a = n.a;
        n.a = temp_a;
    }
}
