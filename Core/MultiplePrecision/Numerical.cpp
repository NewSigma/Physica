/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "CalcBasic.h"
#include <cmath>
#include <iostream>
#include <iomanip>

//Platform dependent
union double_extract {
    double value;
    struct {
        unsigned int low : 32;
        unsigned int high : 20;
        unsigned int exp : 11;
        unsigned int sign : 1;
    } structure;
};
////////////////////////////////Numerical////////////////////////////////
Numerical::Numerical(NumericalUnit* byte, int length, int power, NumericalUnit a) noexcept : byte(byte), length(length), power(power), a(a) {}

Numerical::Numerical(const Numerical& n) noexcept : length(n.length), power(n.power), a(n.a) {
    int size = getSize();
    byte = reinterpret_cast<NumericalUnit*>(malloc(size * sizeof(NumericalUnit)));
    memcpy(byte, n.byte, size * sizeof(NumericalUnit));
}

Numerical::Numerical(Numerical&& n) noexcept : byte(n.byte), length(n.length), power(n.power), a(n.a) {
    n.byte = nullptr;
}

Numerical::Numerical(const Numerical* n) noexcept : Numerical(*n) {}

Numerical::Numerical(SignedNumericalUnit unit, NumericalUnit a) noexcept
: byte(reinterpret_cast<unsigned long*>(malloc(sizeof(NumericalUnit)))), power(0), a(a) {
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
    //10 is the max precision.
    return os << std::setprecision(10) << std::to_string(double(n)) << "\tLength = "
    << n.getSize() << "\tPower = " << n.power << "\tAccuracy = " << (int)n.a << std::setprecision(6); //6 is the default precision.
}

Numerical Numerical::operator<<(int bits) {
    size_t bytes = getSize() * sizeof(NumericalUnit);
    auto new_byte = reinterpret_cast<NumericalUnit*>(malloc(bytes));
    mempcpy(new_byte, byte, bytes);
    return Numerical(new_byte, length, power + bits, a);
}

Numerical Numerical::operator>>(int bits) {
    size_t bytes = getSize() * sizeof(NumericalUnit);
    auto new_byte = reinterpret_cast<NumericalUnit*>(malloc(bytes));
    mempcpy(new_byte, byte, bytes);
    return Numerical(new_byte, length, power - bits, a);
}

NumericalUnit Numerical::operator[](unsigned int index) const {
    return byte[index];
}

Numerical& Numerical::operator= (const Numerical& n) {
    if(this == &n)
        return *this;
    length = n.length;
    int size = getSize();
    free(byte);
    byte = reinterpret_cast<NumericalUnit*>(malloc(size * sizeof(NumericalUnit)));
    memcpy(byte, n.byte, size * sizeof(NumericalUnit));
    power = n.power;
    a = n.a;
    return *this;
}

Numerical& Numerical::operator=(Numerical&& n) noexcept {
    free(byte);
    byte = n.byte;
    n.byte = nullptr;
    length = n.length;
    power = n.power;
    a = n.a;
    return *this;
}
/*
 * Using Newton's method
 */
Numerical Numerical::operator^ (const Numerical& n) const {
    if(isPositive()) {
        if(n.isInteger()) {
            Numerical n2_copy(n);

            Numerical result = getOne();
            if(n.isNegative()) {
                Numerical temp = reciprocal(*this);
                while(n2_copy != basicConst->get_0()) {
                    n2_copy -= basicConst->get_1();
                    result *= temp;
                }
            }
            else {
                while(n2_copy != basicConst->get_0()) {
                    n2_copy -= basicConst->get_1();
                    result *= *this;
                }
            }
            return result;
        }
        else {
            Numerical result = getOne();
            Numerical temp_1 = ln(*this);
            temp_1 *= n;
            temp_1 += basicConst->get_1();
            bool go_on;
            do {
                Numerical temp_2 = ln(result);
                temp_2.toOpposite();
                temp_2 += temp_1;
                temp_2 *= result;
                go_on = temp_2 == result;
                result = std::move(temp_2);
            } while(go_on);
            ++result.a;
            return result;
        }
    }
    else if(isZero() && !n.isZero())
        return getZero();
    else
        qFatal("Can not resolve power of a minus value.");
}

Numerical Numerical::operator-() const {
    size_t bytes = getSize() * sizeof(NumericalUnit);
    auto new_byte = reinterpret_cast<NumericalUnit*>(malloc(bytes));
    mempcpy(new_byte, byte, bytes);
    return Numerical(new_byte, -length, power, a);
}

Numerical operator+ (const Numerical& n1, const Numerical& n2) {
    Numerical result = add(n1, n2);
    cutLength(result);
    if(n1.getA() != 0 || n2.getA() != 0) {
        if(n1.getA() == 0)
            result.applyError(getAccuracy(n2));
        else if(n2.getA() == 0)
            result.applyError(getAccuracy(n1));
        else
            result.applyError(add(getAccuracy(n1), getAccuracy(n2)));
    }
    return result;
}

Numerical operator- (const Numerical& n1, const Numerical& n2) {
    Numerical result = sub(n1, n2);
    cutLength(result);
    if(n1.getA() != 0 || n2.getA() != 0) {
        if(n1.getA() == 0)
            result.applyError(getAccuracy(n2));
        else if(n2.getA() == 0)
            result.applyError(getAccuracy(n1));
        else
            result.applyError(add(getAccuracy(n1), getAccuracy(n2)));
    }
    return result;
}

Numerical operator* (const Numerical& n1, const Numerical& n2) {
    Numerical result = mul(n1, n2);
    cutLength(result);
    if(n1.getA() != 0 || n2.getA() != 0) {
        if(n1.getA() == 0)
            result.applyError(mul(n1, getAccuracy(n2)));
        else if(n2.getA() == 0)
            result.applyError(mul(n2, getAccuracy(n1)));
        else {
            Numerical n1_a = getAccuracy(n1);
            Numerical n2_a = getAccuracy(n2);
            Numerical temp1 = mul(n1_a, n2_a);
            Numerical temp2 = add(mul(n1, n2_a), mul(n2, n1_a));
            result.applyError(add(temp1, temp2));
        }
    }
    return result;
}

Numerical operator/ (const Numerical& n1, const Numerical& n2) {
    Numerical result = div(n1, n2);
    if(n1.getA() != 0 || n2.getA() != 0) {
        if(n1.getA() == 0) {
            Numerical n2_a = getAccuracy(n2);
            Numerical temp_1 = mul(n1, n2_a);
            Numerical temp_2 = mul(n2, sub(n2, n2_a));
            result.applyError(div(temp_1, temp_2));
        }
        else if(n2.getA() == 0)
            result.applyError(div(getAccuracy(n1), n2));
        else {
            Numerical n2_a = getAccuracy(n2);
            Numerical temp_1 = add(mul(n1, n2_a), mul(n2, getAccuracy(n1)));
            Numerical temp_2 = mul(n2, sub(n2, n2_a));
            result.applyError(div(temp_1, temp_2));
        }
    }
    return result;
}

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
    if(n1.getPower() != n2.getPower())
        return false;
    //Here length may not equal n.length because we define numbers like 1.0 and 1.00 are equal.
    if((n1.getLength() ^ n2.getLength()) < 0) // NOLINT(hicpp-signed-bitwise)
        return false;
    if(n1.getA() != n2.getA())
        return false;
    Numerical temp = n1 - n2;
    return temp.isZero();
}

bool operator!= (const Numerical& n1, const Numerical& n2) {
    return !(n1 == n2);
}
//Return accuracy in class Numerical.
Numerical getAccuracy(const Numerical& n) {
    auto b = reinterpret_cast<NumericalUnit*>(malloc(sizeof(NumericalUnit)));
    b[0] = n.getA();
    return Numerical(b, 1, n.getPower() - n.getSize() + 1);
}
//Add error to this and adjust this->length as well as this-> byte.
Numerical& Numerical::applyError(const Numerical& error) {
    if(!error.isZero()) {
        int size = getSize();
        int copy = size;
        int temp = power - size + 1 - error.power;
        NumericalUnit copy_a = a;
        if(temp <= 0) {
            if(temp < 0) {
                auto b = reinterpret_cast<NumericalUnit*>(malloc(sizeof(NumericalUnit)));
                b[0] = a;
                Numerical error_1(b, 1, power - size + 1);
                error_1 = add(error_1, error);
                size += temp;
                a += error_1.byte[error_1.getSize() - 1];
            }
            else
                a += error.byte[error.getSize() - 1];
        }

        if(a < copy_a) {
            a = 1;
            --size;
        }

        if(size < 1) {
            qWarning("Accumulated too many errors.");
            power += 1 - size;
            size = 1;
            if(length < 0)
                size = -size;
            length = size;
            return *this;
        }

        if(size < copy) {
            auto new_byte = reinterpret_cast<NumericalUnit*>(malloc(size * sizeof(NumericalUnit)));
            memcpy(new_byte, byte + copy - size, size * sizeof(NumericalUnit));
            free(byte);
            byte = new_byte;
        }
        if(length < 0)
            size = -size;
        length = size;
    }
    return *this;
}