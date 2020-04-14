/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "Numerical.h"
#include "ElementaryFunction.h"
#include <cstring>
#include <cmath>
#include <QtCore/qlogging.h>
#include <iomanip>
#include "SystemBits.h"

extern const MathConst* mathConst;
/*
 *
 * Useful formulas:
 * (1) n.ed - n.power
 *     Number of digits between the unit and the last digit. (Both included)
 *
 * Operators that do not need to free memory : = += -= *= /= > <
 * Operators that need memory free : new + - * / toNumericalber() randomNumericalber()
 *
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
////////////////////////////////Numerical////////////////////////////////
Numerical::Numerical(unsigned long* byte, int length, int power, unsigned long a) : byte(byte), length(length), power(power), a(a) {}

Numerical::Numerical(const Numerical& n) : length(n.length), power(n.power), a(n.a) {
    int size = getSize();
    byte = (unsigned long*)malloc(size * sizeof(long));
    memcpy(byte, n.byte, size * sizeof(long));
}

Numerical::Numerical(Numerical&& n) noexcept : byte(n.byte), length(n.length), power(n.power), a(n.a) {
    n.byte = nullptr;
}

Numerical::Numerical(const Numerical* n) : Numerical(*n) {}

#pragma clang diagnostic push
#pragma ide diagnostic ignored "hicpp-signed-bitwise"
Numerical::Numerical(double d, unsigned long a) : a(a) {
    if(d == 0) {
        byte = (unsigned long*)malloc(sizeof(long));
        length = 1;
        byte[0] = power = 0;
        return;
    }
    bool negative = d < 0;
    if(negative)
        d = -d;

    power = 0;
    while(d > ULONG_MAX) {
        d /= ULONG_MAX;
        ++power;
    }

    while(d < 1) {
        d *= ULONG_MAX;
        --power;
    }
#ifdef PHYSICA_64BIT
    length = 2;
    byte = (unsigned long*)malloc(length * sizeof(long));
    auto integer = (unsigned long)d;
    byte[1] = integer;
    d -= integer;
    byte[0] = (unsigned long)(d * ULONG_MAX);
#endif

#ifdef PHYSICA_32BIT
    length = 3;
    byte = (unsigned long*)malloc(length * sizeof(long));
    Numerical integer = (unsigned long)d;
    d -= integer;
    byte[2] = integer;
    d *= ULONG_MAX;
    integer = (unsigned long)d;
    d -= integer;
    byte[1] = integer;
    byte[0] = (unsigned long)(d * ULONG_MAX);
#endif
    if(negative)
        length = -length;
}
#pragma clang diagnostic pop
/*
 * May not be very accurate.
 */
Numerical::Numerical(const char* s, unsigned long a) : Numerical(strtod(s, nullptr), a) {}

Numerical::Numerical(const wchar_t* s, unsigned long a) : a(a) {
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

Numerical::Numerical(const std::string& s, unsigned long a) : Numerical(s.c_str(), a) {}

Numerical::Numerical(const std::wstring& s, unsigned long a) : Numerical(s.c_str(), a) {}

Numerical::~Numerical() { free(byte); }
/*
 * Not very accurate.
 */
Numerical::operator double() const {
    int lastIndex = getSize() - 1;
    double d = 0;
    double base = pow(ULONG_MAX, power - lastIndex);
    for(int i = 0; i < lastIndex; ++i) {
        d += base * (double)byte[i];
        base *= ULONG_MAX;
    }
    d += base * (double)byte[lastIndex];

    return d;
}

std::ostream& operator<<(std::ostream& os, const Numerical& n) {
    return os << std::setprecision(53) << (n.length < 0 ? '-' : ' ') << std::to_string(double(n)) << "\tLength = "
    << n.getSize() << "\tPower = " << n.power << "\tAccuracy = " << (int)n.a << std::setprecision(6);
}

void Numerical::operator<<(int bits) {
    power += bits;
}

void Numerical::operator>>(int bits) {
    power -= bits;
}

unsigned long Numerical::operator[](unsigned int index) const {
    return byte[index];
}

Numerical& Numerical::operator= (const Numerical& n) {
    if(this == &n)
        return *this;
    length = n.length;
    int size = getSize();
    free(byte);
    byte = (unsigned long*)malloc(size * sizeof(long));
    memcpy(byte, n.byte, size * sizeof(long));
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

Numerical Numerical::operator+ (const Numerical& n) const {
    Numerical result = add(*this, n);
    cutLength(result);
    if(a != 0 || n.a != 0) {
        if(a == 0)
            result.applyError(n.getAccuracy());
        else if(n.a == 0)
            result.applyError(getAccuracy());
        else
            result.applyError(add(getAccuracy(), n.getAccuracy()));
    }
    return result;
}

Numerical Numerical::operator- (const Numerical& n) const {
    Numerical result = sub(*this, n);
    cutLength(result);
    if(a != 0 || n.a != 0) {
        if(a == 0)
            result.applyError(n.getAccuracy());
        else if(n.a == 0)
            result.applyError(getAccuracy());
        else
            result.applyError(add(getAccuracy(), n.getAccuracy()));
    }
    return result;
}

Numerical Numerical::operator* (const Numerical& n) const {
    Numerical result = mul(*this, n);
    cutLength(result);
    if(a != 0 || n.a != 0) {
        if(a == 0)
            result.applyError(mul(*this, n.getAccuracy()));
        else if(n.a == 0)
            result.applyError(mul(n, getAccuracy()));
        else {
            Numerical n1_a = getAccuracy();
            Numerical n2_a = n.getAccuracy();
            Numerical temp1 = mul(n1_a, n2_a);
            Numerical temp2 = add(mul(*this, n2_a), mul(n, n1_a));
            result.applyError(add(temp1, temp2));
        }
    }
    return result;
}

Numerical Numerical::operator/ (const Numerical& n) const {
    Numerical result = div(*this, n);
    if(a != 0 || n.a != 0) {
        if(a == 0) {
            Numerical n2_a = n.getAccuracy();
            Numerical temp_1 = mul(*this, n2_a);
            Numerical temp_2 = mul(n, sub(n, n2_a));
            result.applyError(div(temp_1, temp_2));
        }
        else if(n.a == 0)
            result.applyError(div(getAccuracy(), n));
        else {
            Numerical n2_a = n.getAccuracy();
            Numerical temp_1 = add(mul(*this, n2_a), mul(n, getAccuracy()));
            Numerical temp_2 = mul(n, sub(n, n2_a));
            result.applyError(div(temp_1, temp_2));
        }
    }
    return result;
}
/*
 * Using Newton's method
 * May return nullptr.
 */
Numerical Numerical::operator^ (const Numerical& n) const {
    if(Q_UNLIKELY(isZero())) {
        if(!n.isZero())
            return getZero();
    }
    else if(isPositive()) {
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
            result.a += 1;
            return result;
        }
    }
}

bool Numerical::operator> (const Numerical& n) const {
    //Judge from sign.
    if(isPositive()) {
        if(!n.isPositive())
            return true;
    }
    else if(isZero())
        return n.isNegative();
    else {
        if(!n.isNegative())
            return false;
    }
    //If we cannot get a result, judge from power
    bool result;
    if(power > n.power)
        result = true;
    else if(power < n.power)
        result = false;
    else {
        //The only method left.
        Numerical subtract = *this - n;
        result = subtract.isPositive();
    }
    return result;
}

bool Numerical::operator< (const Numerical& n) const {
    //Judge from sign.
    if(isPositive()) {
        if(!n.isPositive())
            return false;
    }
    else if(isZero())
        return n.isPositive();
    else {
        if(!n.isNegative())
            return true;
    }
    //If we cannot get a result, judge from power
    bool result;
    if(power > n.power)
        result = false;
    else if(power < n.power)
        result = true;
    else {
        //The only method left.
        Numerical subtract = *this - n;
        result = subtract.isNegative();
    }
    return result;
}

bool Numerical::operator>= (const Numerical& n) const {
    return !(*this < n);
}

bool Numerical::operator<= (const Numerical& n) const {
    return !(*this > n);
}

bool Numerical::operator== (const Numerical& n) const {
    if(power != n.power)
        return false;
    //Here length may not equal n.length because we define numbers like 1.0 and 1.00 are equal.
    if((length ^ n.length) < 0) // NOLINT(hicpp-signed-bitwise)
        return false;
    if(a != n.a)
        return false;
    Numerical temp = *this - n;
    return temp.isZero();
}

bool Numerical::operator!= (const Numerical& n) const {
    return !(*this == n);
}

Numerical Numerical::operator-() const {
    size_t bytes = getSize() * sizeof(long);
    auto new_byte = (unsigned long*)malloc(bytes);
    mempcpy(new_byte, byte, bytes);
    return Numerical(new_byte, length, power, a);
};
//Return accuracy in class Numerical.
Numerical Numerical::getAccuracy() const {
    auto b = (unsigned long*)malloc(sizeof(long));
    b[0] = a;
    return Numerical(b, 1, power - getSize() + 1);
}
//Add error to this and adjust this->length as well as this-> byte.
Numerical& Numerical::applyError(const Numerical& error) {
    if(!error.isZero()) {
        int size = getSize();
        int copy = size;
        int temp = power - size + 1 - error.power;
        unsigned long copy_a = a;
        if(temp <= 0) {
            if(temp < 0) {
                auto b = (unsigned long*)malloc(sizeof(long));
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
            auto new_byte = (unsigned long*)malloc(size * sizeof(long));
            memcpy(new_byte, byte + copy - size, size * sizeof(long));
            free(byte);
            byte = new_byte;
        }
        if(length < 0)
            size = -size;
        length = size;
    }
    return *this;
}

void printElements(const Numerical& n) {
    int size = n.getSize();
    for(int i = 0; i < size; ++i)
        std::cout << n[i] << ' ';
}