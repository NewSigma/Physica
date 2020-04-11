/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "Numerical.h"
#include "ElementaryFunction.h"
#include <cstring>
#include <cmath>
#include <QtCore/qlogging.h>
#include <vector>
#include <iomanip>
#include "CalcBasic.h"
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
Numerical::Numerical(unsigned long* byte, signed char length, int power, unsigned char a) : byte(byte), power(power), length(length), a(a) {}

Numerical::Numerical(const Numerical& n) : power(n.power), length(n.length), a(n.a) {
    auto size = getSize();
    byte = (unsigned long*)malloc(size * sizeof(long));
    memcpy(byte, n.byte, size * sizeof(long));
}

Numerical::Numerical(const Numerical* n) : Numerical(*n) {}

#pragma clang diagnostic push
#pragma ide diagnostic ignored "hicpp-signed-bitwise"
Numerical::Numerical(double d, unsigned char acc) : a(acc) {
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
    auto integer = (unsigned long)d;
    d -= integer;
    byte[2] = integer;
    d *= ULONG_MAX;
    integer = (unsigned long)d;
    d -= integer;
    byte[1] = integer;
    byte[0] = (unsigned long)(d * ULONG_MAX);
#endif
    if(negative)
        length = (signed char)-length;
}
#pragma clang diagnostic pop
/*
 * May not be very accurate.
 */
Numerical::Numerical(const char* s, unsigned char a) : Numerical(strtod(s, nullptr), a) {}

Numerical::Numerical(const wchar_t* s, unsigned char a) : a(a) {
    signed char size = wcslen(s);
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

Numerical::Numerical(const std::string& s, unsigned char a) : Numerical(s.c_str(), a) {}

Numerical::Numerical(const std::wstring& s, unsigned char a) : Numerical(s.c_str(), a) {}

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
//Move operator, which moves this to n.
void Numerical::operator<<(Numerical& n) {
    free(byte);
    byte = n.byte;
    length = n.length;
    power = n.power;
    a = n.a;
    n.byte = nullptr;
    delete &n;
}

Numerical& Numerical::operator= (const Numerical& n) {
    if(this == &n)
        return *this;
    length = n.length;
    int size = getSize();
    byte = (unsigned long*)realloc(byte, size * sizeof(long));
    memcpy(byte, n.byte, size * sizeof(long));
    power = n.power;
    a = n.a;
    return *this;
}

Numerical* Numerical::operator+ (const Numerical& n) const {
    auto result = add(*this, n);
    cutLength(result);
    if(a != 0 || n.a != 0) {
        Numerical* error;
        if(a == 0)
            error = n.getAccuracy();
        else if(n.a == 0)
            error = getAccuracy();
        else {
            auto n1_a = getAccuracy();
            auto n2_a = n.getAccuracy();
            error = add(*n1_a, *n2_a);
            delete n1_a;
            delete n2_a;
        }
        result->applyError(error);
        delete error;
    }
    return result;
}

Numerical* Numerical::operator- (const Numerical& n) const {
    auto result = subtract(*this, n);
    cutLength(result);
    if(a != 0 || n.a != 0) {
        Numerical* error;
        if(a == 0)
            error = n.getAccuracy();
        else if(n.a == 0)
            error = getAccuracy();
        else {
            auto n1_a = getAccuracy();
            auto n2_a = n.getAccuracy();
            error = add(*n1_a, *n2_a);
            delete n1_a;
            delete n2_a;
        }
        result->applyError(error);
        delete error;
    }
    return result;
}

Numerical* Numerical::operator* (const Numerical& n) const {
    auto result = multiply(*this, n);
    cutLength(result);
    if(a != 0 || n.a != 0) {
        Numerical* error;
        if(a == 0) {
            auto n2_a = n.getAccuracy();
            error = multiply(*this, *n2_a);
            delete n2_a;
        }
        else if(n.a == 0) {
            auto n1_a = getAccuracy();
            error = multiply(n, *n1_a);
            delete n1_a;
        }
        else {
            auto n1_a = getAccuracy();
            auto n2_a = n.getAccuracy();
            auto temp_1 = multiply(*this, *n2_a);
            auto temp_2 = multiply(n, *n1_a);
            auto temp_3 = multiply(*n1_a, *n2_a);
            auto temp_4 = add(*temp_1, *temp_2);
            error = add(*temp_3, *temp_4);
            delete n1_a;
            delete n2_a;
            delete temp_1;
            delete temp_2;
            delete temp_3;
            delete temp_4;
        }
        result->applyError(error);
        delete error;
    }
    return result;
}

Numerical* Numerical::operator/ (const Numerical& n) const {
    auto result = divide(*this, n);
    if(result != nullptr) {
        if(a != 0 || n.a != 0) {
            Numerical* error;
            if(a == 0) {
                auto n2_a = n.getAccuracy();
                auto temp_1 = multiply(*this, *n2_a);
                auto temp_2 = subtract(n, *n2_a);
                auto temp_3 = multiply(n, *temp_2);
                error = divide(*temp_1, *temp_3);
                delete n2_a;
                delete temp_1;
                delete temp_2;
                delete temp_3;
            }
            else if(n.a == 0) {
                auto n1_a = getAccuracy();
                error = divide(*n1_a, n);
                delete n1_a;
            }
            else {
                auto n1_a = getAccuracy();
                auto n2_a = n.getAccuracy();
                auto temp_1 = multiply(*this, *n2_a);
                auto temp_2 = multiply(n, *n1_a);
                auto temp_3 = add(*temp_1, *temp_2);
                auto temp_4 = subtract(n, *n2_a);
                auto temp_5 = multiply(n, *temp_4);
                error = divide(*temp_3, *temp_5);
                delete n1_a;
                delete n2_a;
                delete temp_1;
                delete temp_2;
                delete temp_3;
                delete temp_4;
                delete temp_5;
            }
            result->applyError(error);
            delete error;
        }
    }
    return result;
}
/*
 * Using Newton's method
 * May return nullptr.
 */
Numerical* Numerical::operator^ (const Numerical& n) const {
    Numerical* result = nullptr;
    if(Q_UNLIKELY(isZero())) {
        if(!n.isZero())
            result = getZero();
    }
    else if(isPositive()) {
        if(n.isInteger()) {
            auto n2_copy = new Numerical(n);

            result = getOne();
            if(n.isNegative()) {
                auto temp = reciprocal(*this);
                while(*n2_copy != basicConst->get_0()) {
                    *n2_copy -= basicConst->get_1();
                    *result *= *temp;
                }
                delete temp;
            }
            else {
                while(*n2_copy != basicConst->get_0()) {
                    *n2_copy -= basicConst->get_1();
                    *result *= *this;
                }
            }
            delete n2_copy;
        }
        else {
            auto temp_result = getOne();
            auto temp_1 = ln(*this);
            *temp_1 *= n;
            *temp_1 += basicConst->get_1();
            bool go_on;
            do {
                result = ln(*temp_result);
                result->length = (signed char)-result->length;
                *result += *temp_1;
                *result *= *temp_result;
                go_on = *result == *temp_result;
                delete temp_result;
                temp_result = result;
            } while(go_on);
            result->a += 1;
            delete temp_1;
        }
    }
    return result;
}
/*
 * Warning: *this can not be temp object.
 */
void Numerical::operator+= (const Numerical& n) {
    Numerical* p_result = *this + n;
    free(byte);
    byte = p_result->byte;
    length = p_result->length;
    power = p_result->power;
    a = p_result->a;
    p_result->byte = nullptr;
    delete p_result;
}

void Numerical::operator-= (const Numerical& n) {
    Numerical* p_result = *this - n;
    free(byte);
    byte = p_result->byte;
    length = p_result->length;
    power = p_result->power;
    a = p_result->a;
    p_result->byte = nullptr;
    delete p_result;
}

void Numerical::operator*= (const Numerical& n) {
    Numerical* p_result = *this * n;
    free(byte);
    byte = p_result->byte;
    length = p_result->length;
    power = p_result->power;
    a = p_result->a;
    p_result->byte = nullptr;
    delete p_result;
}
/*
 * n2 mustn't be zero.
 */
void Numerical::operator/= (const Numerical& n) {
    Numerical* p_result = *this / n;
    free(byte);
    byte = p_result->byte;
    length = p_result->length;
    power = p_result->power;
    a = p_result->a;
    p_result->byte = nullptr;
    delete p_result;
}

void Numerical::operator^= (const Numerical& n) {
    Numerical* p_result = *this ^ n;
    free(byte);
    byte = p_result->byte;
    length = p_result->length;
    power = p_result->power;
    a = p_result->a;
    p_result->byte = nullptr;
    delete p_result;
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
        Numerical* p_result = *this - n;
        result = p_result->isPositive();
        delete p_result;
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
        Numerical* p_result = *this - n;
        result = p_result->isNegative();
        delete p_result;
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
    auto temp = *this - n;
    if(!temp->isZero())
        return false;
    delete temp;
    return true;
}

bool Numerical::operator!= (const Numerical& n) const {
    return !(*this == n);
}

Numerical* Numerical::operator- () const {
    auto result = new Numerical(this);
    result->length = (signed char)-result->length;
    return result;
}
//Return accuracy in class Numerical.
Numerical* Numerical::getAccuracy() const {
    auto b = (unsigned long*)malloc(sizeof(long));
    b[0] = a;
    return new Numerical(b, 1, power - getSize() + 1);
}
//Return this + accuracy
Numerical* Numerical::getMaximum() const {
    auto acc = getAccuracy();
    auto result = add(*this, *acc);
    delete acc;
    return result;
}
//Return this - accuracy
Numerical* Numerical::getMinimum() const {
    auto acc = getAccuracy();
    auto result = subtract(*this, *acc);
    delete acc;
    return result;
}
//Add error to this and adjust this->length as well as this-> byte.
void Numerical::applyError(const Numerical* error) {
    if(!error->isZero()) {
        int size = getSize();
        int copy = size;
        int temp = power - size + 1 - error->power;
        unsigned char copy_a = a;
        if(temp <= 0) {
            if(temp < 0) {
                auto b = (unsigned long*)malloc(sizeof(long));
                b[0] = a;
                auto error_1 = Numerical(b, 1, power - size + 1);
                error_1 << *add(error_1, *error);
                size += temp;
                a += error_1.byte[error_1.getSize() - 1];
            }
            else
                a += error->byte[error->getSize() - 1];
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
            length = (signed char)size;
            return;
        }

        if(size < copy) {
            auto new_byte = (unsigned long*)malloc(size * sizeof(long));
            memcpy(new_byte, byte + copy - size, size * sizeof(long));
            free(byte);
            byte = new_byte;
        }
        if(length < 0)
            size = -size;
        length = (signed char)size;
    }
}

void Numerical::printElements() const {
    for(int i = 0; i < getSize(); ++i)
        std::cout << byte[i] << ' ';
}