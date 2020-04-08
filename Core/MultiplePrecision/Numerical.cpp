/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "Numerical.h"
#include "Solve.h"
#include <cstring>
#include <cmath>
#include <QtCore/qlogging.h>

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
Numerical::Numerical() : Numerical(nullptr, 0, 0) {}

Numerical::Numerical(const Numerical& n) : power(n.power), length(n.length), a(n.a) {
    auto size = getSize();
    byte = (unsigned char*)malloc(size * sizeof(char));
    memcpy(byte, n.byte, size * sizeof(char));
}

Numerical::Numerical(unsigned char* byte, signed char length, int power, unsigned char a) : byte(byte), power(power), length(length), a(a) {}

Numerical::Numerical(const Numerical* n) : Numerical(*n) {}
/*
 * May not be very accurate.
 */
Numerical::Numerical(double d, unsigned char acc) {
    bool negative = d < 0;
    d = negative ? -d : d;
    auto copy_d = d;

    power = length = 0;
    if(d < 1) {
        do {
            d *= 10;
            --power;
        } while(d < 1);
    }
    else if(d >= 10) {
        do {
            d /= 10;
            ++power;
        } while(d >= 10);
    }

    do {
        copy_d -= double(int(copy_d));
        copy_d *= 10;
        ++length;
    } while(copy_d != 0 && length <= basicConst->getGlobalPrecision());

    byte = (unsigned char*)malloc(length * sizeof(char));
    int copy = length;
    while(--copy) {
        byte[copy] = (char)d;
        d -= double(int(d));
        d *= 10;
    }
    byte[0] = (char)d;

    if(negative)
        length = (signed char)-length;
    a = acc;
}

Numerical::Numerical(const char* s, unsigned char acc) {
    signed char size = strlen(s);
    byte = new unsigned char[size];
    bool negative = false;
    a = acc;
    signed char index = 0, id_byte = 0, start_eff_num = 0, point_id = size;
    for(; index < size; ++index) {
        switch(s[index]) {
            case '-':
                negative = true;
                break;
            case '.':
                point_id = index;
            case '0':
                break;
            case '1' ... '9':
                start_eff_num = index;
                goto double_break;
            default:
                qCritical("Failed to initialize Numerical.");
                byte = nullptr;
                power = a = length = 0;
                return;
        }
    }
    double_break:
    for(; index < size; ++index) {
        if(s[index] != '.') {
            byte[id_byte] = s[index] - '0';
            ++id_byte;
        }
        else {
            point_id = index;
        }
    }
    length = id_byte > basicConst->getGlobalPrecision() ? basicConst->getGlobalPrecision() : id_byte;
    power = point_id - start_eff_num + (point_id > start_eff_num ? -1 : 0);
    byte = (unsigned char*)realloc(byte, length);
    reverse(byte, getSize());
    if(negative)
        length = (signed char)-length;
}

Numerical::Numerical(const wchar_t* s, unsigned char acc) {
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
    a = acc;
}

Numerical::Numerical(const std::string& s, unsigned char acc) : Numerical(s.c_str(), acc) {}

Numerical::Numerical(const std::wstring& s, unsigned char acc) : Numerical(s.c_str(), acc) {}

Numerical::~Numerical() { free(byte); }

std::string Numerical::toString() const {
    std::string result;
    int size = getSize();
    if(byte != nullptr) {
        if(isNegative())
            result.push_back('-');
        if(abs(power) > basicConst->getGlobalPrecision()) {
            result.push_back(byte[size - 1] + '0');
            --size;
            result.push_back('.');
            while(size--)
                result.push_back(byte[size] + '0');
            result += "×10^";
            result += std::to_string(power);
        }
        else if(power < 0) {
            result += "0.";
            for(int i = -power - 1; i > 0; --i)
                result.push_back('0');
            while(size--)
                result.push_back(byte[size] + '0');
        }
        else {
            if(size > power) {
                for (int i = 1; i <= power + 1; ++i)
                    result.push_back(byte[size - i] + '0');
                result.push_back('.');
                for(int i = power + 2; i <= size; ++i)
                    result.push_back(byte[size - i] + '0');
            }
            else {
                while(size--)
                    result.push_back(byte[size] + '0');
                for(int i = power - size + 1; i > 0; --i)
                    result.push_back('0');
            }
        }
    }
    return result;
}
/*
 * Not very accurate either.
 */
Numerical::operator double() const {
    double result_integer = 0;
    double result_float = 0;
    int size = getSize();

    int temp_index = power + 1;
    for(int i = 0; i < temp_index; ++i) {
        result_integer *= 10;
        result_integer += byte[i];
    }

    while(temp_index < size) {
        result_float += byte[size - temp_index - 1];
        result_float *= 10;
        ++temp_index;
    }
    result_float /= pow(10,size - power);

    return result_integer + result_float;
}

std::ostream& operator<<(std::ostream& os, const Numerical& n) {
    int size = n.getSize();
    os << n.toString() << "\tLength = " << size << "\tPower = " << n.power << "\tAccuracy = ";

    int temp = n.power - size +1;
    if(n.a == 0)
        os << "0.";
    else {
        if(temp < 0) {
            os << "0.";
            for(; temp < -1; ++temp)
                os << '0';
            os << (int)n.a;
        }
        else {
            os << (int)n.a;
            for(; temp > 0; --temp)
                os << '0';
            os << '.';
        }
    }
    return os;
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
    byte = (unsigned char*)realloc(byte, size * sizeof(char));
    memcpy(byte, n.byte, size * sizeof(char));
    power = n.power;
    a = n.a;
    return *this;
}

Numerical* Numerical::operator+ (const Numerical& n) const {
    auto result = add(*this, n);
    result->a = cutLength(result);
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
    result->a = cutLength(result);
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
    result->a = cutLength(result);
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
            auto boolean = result->applyError(error);
            delete error;
            if(boolean)
                return nullptr;
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
    *this << *p_result;
    delete p_result;
}

bool Numerical::operator> (const Numerical& n) const {
    //Judge from sign.
    if(isPositive()) {
        if(n.isZero() || n.isNegative())
            return true;
    }
    else if(isZero())
        return n.isNegative();
    else {
        if(n.isPositive() || n.isZero())
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
        if(n.isZero() || n.isNegative())
            return false;
    }
    else if(isZero())
        return n.isPositive();
    else {
        if(n.isPositive() || n.isZero())
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
    auto b = (unsigned char*)malloc(sizeof(char));
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
bool Numerical::applyError(const Numerical* error) {
    int size = getSize();
    int copy = size;
    int temp = power - size + 1 - error->power;
    if(temp <= 0) {
        if(temp < 0) {
            auto b = (unsigned char*)malloc(sizeof(char));
            b[0] = a;
            auto error_1 = new Numerical(b, 1, power - size + 1);
            *error_1 << *add(*error_1, *error);
            size += temp;
            a += error_1->byte[error_1->getSize() - 1];
            delete error_1;
        }
        else
            a += error->byte[error->getSize() - 1];
    }

    if(a > 9) {
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
        return true;
    }

    if(size < copy) {
        auto new_byte = (unsigned char*)malloc(size * sizeof(char));
        memcpy(new_byte, byte + copy - size, size * sizeof(char));
        free(byte);
        byte = new_byte;
    }
    if(length < 0)
        size = -size;
    length = (signed char)size;
    return false;
}
////////////////////////////////Helper functions/////////////////////////////////////
//reverse the order of elements of arr
void reverse(unsigned char* arr, int length) {
    unsigned char temp;
    for(int i = 0; i < length / 2; ++i) {
        temp = arr[i];
        arr[i] = arr[length - i - 1];
        arr[length - i - 1] = temp;
    }
}
//////////////////////////////Process functions////////////////////////////////////////
/*
 * The following four functions simply calculate the result while operator functions will
 * consider the accuracy.
 */
Numerical* add (const Numerical& n1, const Numerical& n2) {
    Numerical* result;
    if(n1.isZero())
        result = new Numerical(n2);
    else if(n2.isZero())
        result = new Numerical(n1);
    else if ((n1.length ^ n2.length) < 0) { // NOLINT(hicpp-signed-bitwise)
        Numerical* shallow_copy;
        if (n1.length > 0) {
            shallow_copy = new Numerical(n2.byte, n2.length, n2.power);
            result = n1 - *shallow_copy;
            shallow_copy->byte = nullptr;
        }
        else {
            shallow_copy = new Numerical(n1.byte, n1.length, n1.power);
            result = n2 - *shallow_copy;
            shallow_copy->byte = nullptr;
        }
        delete shallow_copy;
    }
    else {
        const Numerical* big;
        const Numerical* small;
        if (n1.power > n2.power) {
            big = &n1;
            small = &n2;
        }
        else {
            big = &n2;
            small = &n1;
        }
        int bigSize = big->getSize();
        int smallSize = small->getSize();
        //Estimate the ed of result first, will calculate it accurately later.
        signed char length = (signed char)(big->power + std::max(bigSize - big->power, smallSize - small->power));
        auto byte = (unsigned char*)malloc(length * sizeof(char));
        memcpy(byte + length - bigSize, big->byte, bigSize * sizeof(char));
        memset(byte, 0, (length - bigSize) * sizeof(char));
        //Add small to big
        int index = length - big->power + small->power - smallSize;
        int lastIndex = smallSize - 1;
        for(int i = 0; i < lastIndex; ++i) {
            byte[index] += small->byte[i];
            //Carry bit;
            if(byte[index] > 9) {
                ++byte[index + 1];
                byte[index] -= 10;
            }
            ++index;
        }
        byte[index] += small->byte[lastIndex];
        //Carry bit
        while(byte[index] > 9 && index != length - 1) {
            byte[index] -= 10;
            ++index;
            ++byte[index];
        }
        ///////////////////////////////////////Get byte, length and power//////////////////////////
        int power = big->power;
        if(byte[0] > 9) {
            ++length;
            ++power;
            byte = (unsigned char*)realloc(byte, length * sizeof(char));
            byte[length - 1] = 1;
            byte[length - 2] -= 10;
        }
        ////////////////////////////////////Out put////////////////////////////////////////
        if(big->length < 0)
            length = (signed char)-length;
        result = new Numerical(byte, length, power);
    }
    return result;
}

Numerical* subtract (const Numerical& n1, const Numerical& n2) {
    Numerical* result;
    ///////////////////////////////////Deal with Binderies////////////////////////////////
    Numerical* shallow_copy = nullptr;
    if(n1.isZero())
        result = -n2;
    else if(n2.isZero())
        result = new Numerical(n1);
    else if (n1.length > 0) {
        if (n2.length < 0) {
            shallow_copy = new Numerical(n2.byte, n2.length, n2.power);
            result = n1 + *shallow_copy;
            shallow_copy->byte = nullptr;
        }
        else {
            const Numerical* big;
            const Numerical* small;
            bool changeSign = false;
            if (n1.power > n2.power) {
                big = &n1;
                small = &n2;
            }
            else {
                changeSign = true;
                big = &n2;
                small = &n1;
            }
            int bigSize = big->getSize();
            int smallSize = small->getSize();
            //Estimate the ed of result first, will calculate it accurately later.
            signed char length = (signed char)(big->power + std::max(bigSize - big->power, smallSize - small->power));
            auto byte = (signed char*)malloc(length * sizeof(char));
            memcpy(byte + length - bigSize, big->byte, bigSize * sizeof(char));
            memset(byte, 0, (length - bigSize) * sizeof(char));

            int index = length - big->power + small->power - smallSize;
            int lastIndex = smallSize - 1;
            //Add small to big
            for(int i = 0; i < lastIndex; ++i) {
                byte[index] = (signed char)(byte[index] - small->byte[i]);
                //Carry bit;
                if(byte[index] < 0) {
                    --byte[index + 1];
                    byte[index] += 10;
                }
                ++index;
            }
            byte[index] = (signed char)(byte[index] - small->byte[lastIndex]);
            //Carry bit
            while(byte[index] < 0 && index != length - 1) {
                byte[index] += 10;
                ++index;
                --byte[index];
            }

            if(byte[length - 1] >= 0) {
                if(changeSign)
                    length = (signed char)-length;
                result = new Numerical((unsigned char*)byte, length, big->power);
            }
            else {
                //If n1 - n2 < 0, we have to change our method.
                free(byte);
                result = subtract(n2, n1);
                result->length = (signed char)-result->length;
            }
            cutZero(result);
        }
    }
    else {
        shallow_copy = new Numerical(n1.byte, n1.length, n1.power);
        if (n2.length > 0) {
            result = *shallow_copy + n2;
            result->length = (signed char)-result->length;
        }
        else {
            auto shallow_copy_1 = new Numerical(n2.byte, n2.length, n2.power);
            result = *shallow_copy_1 - *shallow_copy;
            shallow_copy_1->byte = nullptr;
            delete shallow_copy_1;
        }
        shallow_copy->byte = nullptr;
    }
    ////////////////////////////////////Out put////////////////////////////////////////
    delete shallow_copy;
    return result;
}

Numerical* multiply (const Numerical& n1, const Numerical& n2) {
    Numerical* result;
    if(n1 == basicConst->get_1())
        result = new Numerical(n2);
    else if(n2 == basicConst->get_1())
        result = new Numerical(n1);
    else {
        int size1 = n1.getSize();
        int size2 = n2.getSize();
        int last1 = size1 - 1;
        int last2 = size2 - 1;
        //Estimate the ed of result first. we will calculate it accurately later.
        auto length = (signed char)(size1 + size2 - 1);
        auto lastIndex = length - 1;
        auto byte = (unsigned char*)calloc(length, sizeof(char));
        //i = 1 ... size - 2
        for (int i = 0; i < last1; ++i) {
            for(int j = 0; j < size2; ++j) {
                int index = i + j;
                byte[index] += n1.byte[i] * n2.byte[j];
                //Carry bit.
                if (byte[index] > 9) {
                    int tens = byte[index] / 10;
                    byte[index] -= tens * 10;
                    byte[index + 1] += tens;
                }
            }
        }
        //i = size - 1
        for(int j = 0; j < last2; ++j) {
            int index = size1 - 1 + j;
            byte[index] += n1.byte[last1] * n2.byte[j];
            //Carry bit.
            if (byte[index] > 9) {
                int tens = byte[index] / 10;
                byte[index] -= tens * 10;
                byte[index + 1] += tens;
            }
        }
        byte[last1 + last2] += n1.byte[last1] * n2.byte[last2];
        ///////////////////////////////////////Get byte, length and power//////////////////////////;
        int power = n1.power + n2.power;
        if (byte[lastIndex] > 9) {
            ++length;
            ++power;
            int tens = byte[lastIndex] / 10;
            byte[lastIndex] -= tens * 10;
            byte = (unsigned char*)realloc(byte, length * sizeof(char));
            byte[length - 1] = tens;
        }
        ////////////////////////////////////Out put////////////////////////////////////////
        if((n1.length ^ n2.length) < 0) // NOLINT(hicpp-signed-bitwise)
            length = (signed char)-length;
        result = new Numerical(byte, length, power);
    }
    return result;
}

Numerical* divide (const Numerical& n1, const Numerical& n2) {
    if(Q_UNLIKELY(n2.isZero()))
        qFatal("Divide by zero!");

    Numerical* result;
    if(Q_UNLIKELY(n1.isZero() || n2.isZero()))
        result = getZero();
    else if(Q_UNLIKELY(n2 == basicConst->get_1()))
        result = new Numerical(n1);
    else {
        auto n1_copy = new Numerical(n1);
        auto n2_copy = new Numerical(n2);
        n1_copy->length = (signed char)n1_copy->getSize();
        n2_copy->length = (signed char)n2_copy->getSize();
        ////////////////////////////////Calculate cursory first//////////////////////////////////////
        //Estimate the length of result.
        signed char length = basicConst->getGlobalPrecision();
        //Change n1_copy's power larger than n2_copy, power of the result will change correspondingly.
        int power = n1_copy->power - n2_copy->power - 1;
        n1_copy->power = n2_copy->power + 1;
        auto byte = (unsigned char*)calloc(length, sizeof(char));

        auto n1_copy_old = n1_copy;
        for (int i = length - 1; i >= 0; --i) {
            unsigned char unit = 0;
            while(true) {
                n1_copy = subtract(*n1_copy, *n2_copy);
                if(Q_UNLIKELY(n1_copy->isNegative())) {
                    delete n1_copy;
                    n1_copy = n1_copy_old;
                    break;
                }
                else {
                    ++unit;
                    delete n1_copy_old;
                    n1_copy_old = n1_copy;
                    if(n1_copy->isZero()) {
                        byte[i] = unit;
                        goto double_break;
                    }
                }
            }
            ++n1_copy->power;
            byte[i] = unit;
        }
        double_break:
        delete n1_copy;
        delete n2_copy;
        ////////////////////////////////////Out put////////////////////////////////////////
        if(byte[length - 1] > 9) {
            int tens = byte[length - 1] / 10;
            byte[length - 1] -= tens * 10;
            ++power;
            ++length;
            byte = (unsigned char*)realloc(byte, length * sizeof(char));
            byte[length - 1] = tens;
        }
        if((n1.length ^ n2.length) < 0) // NOLINT(hicpp-signed-bitwise)
            length = (signed char)-length;
        //1 comes from the algorithm
        result = new Numerical(byte, length, power, 1);
    }
    return result;
}
/*
 * If the length of new array is larger than GlobalPrecision, it will be set to GlobalPrecision.
 * Return true if array is cut.
 */
bool cutLength(Numerical* n) {
    bool result = false;
    int size = n->getSize();

    if(size > basicConst->getGlobalPrecision()) {
        result = true;
        int cutFrom = size - basicConst->getGlobalPrecision();
        auto new_byte = (unsigned char*)malloc(basicConst->getGlobalPrecision() * sizeof(char));
        memcpy(new_byte, n->byte + cutFrom, basicConst->getGlobalPrecision() * sizeof(char));
        free(n->byte);
        n->byte = new_byte;
        auto length = basicConst->getGlobalPrecision();
        if(n->length < 0)
            length = -length;
        n->length = length;
    }
    return result;
}
/*
 * Cut zeros from the beginning.
 */
void cutZero(Numerical* n) {
    int size = n->getSize();
    int id = size - 1;
    while(n->byte[id] == 0 && id > 0)
        --id;
    ++id;

    if(id != size) {
        int shorten = size - id;
        n->byte = (unsigned char*)realloc(n->byte, id * sizeof(char));
        size = id;
        if(n->length < 0)
            size = -size;
        n->length = (signed char)size;

        if(n->byte[id - 1] != 0)
            n->power -= shorten;
        else
            n->power = 0;
    }
}
//////////////////////////////Basic Operations////////////////////////////////////////
//Return a real number between 0 and 1.
Numerical* randomNumerical() {
    srand(clock());
    srand(random());
    auto result = new Numerical((double)random());
    *result /= basicConst->getR_MAX();

    return result;
}
//Return a real number lowerBound and upperBound.
Numerical* randomNumerical(Numerical* lowerBound, Numerical* upperBound) {
    Numerical* random = randomNumerical();
    auto result = *lowerBound - *upperBound;
    *random *= *random;
    *result += *lowerBound;
    delete random;

    return result;
}
////////////////////////////////////////Elementary Functions////////////////////////////////////////////
Numerical* reciprocal(const Numerical& n) {
    return basicConst->get_1() / n;
}
/*
 * *_light functions do not consider the error caused by a. For example, sqrt_light does not calculate
 * (sqrt(n + a) - sqrt(n)) for error.
 */
Numerical* sqrt_light(const Numerical& n) {
    if(n.length < 0)
        return nullptr;
    auto MachinePrecision = basicConst->getGlobalPrecision();
    auto copy_n = new Numerical(n);
    //Let n < 1 so as to control error.
    char add_power = 0;
    if(copy_n->power > 0) {
        if(copy_n->power % 2 == 0) {
            add_power = char(copy_n->power / 2 + 1);
            copy_n->power = -2;
        }
        else {
            add_power = char((copy_n->power + 1) / 2);
            copy_n->power = -1;
        }
    }

    Numerical* result = getOne();
    Numerical* temp;
    //3.33 is the big approximate value of ln(10)/ln(2)
    for(int i = 0; i < 3.33 * MachinePrecision; ++i) {
        temp = divide(*copy_n, *result);
        *result += *temp;
        *result /= basicConst->get_2();
        delete temp;
    }
    delete copy_n;
    result->power += add_power;
    result->a = 1;

    return result;
}

Numerical* sqrt(const Numerical& n) {
    auto result  = sqrt_light(n);
    if(n.a != 0) {
        auto n_error = n.getMinimum();
        auto error = sqrt_light(*n_error);
        if(error == nullptr) {
            delete n_error;
            n_error = n.getMaximum();
            error = sqrt_light(*n_error);
        }
        *error -= *result;
        error->length = (signed char)error->getSize();

        result->applyError(error);

        delete n_error;
        delete error;
    }
    return result;
}
//TODO not completed: Use gamma function.
Numerical* factorial(const Numerical& n) {
    if(n.length < 0) {
        qCritical("Cannot solve the factorial of a minus value.");
        return nullptr;
    }

    Numerical* result;
    if(n.isInteger()) {
        result = getOne();
        auto temp = getOne();
        while(*temp < n) {
            *temp += basicConst->get_1();
            *result *= *temp;
        }
        delete temp;
        return result;
    }
    return nullptr;
}

Numerical* ln_light(const Numerical& n) {
    if(!n.isPositive())
        return nullptr;
    auto result = getZero();
    if(n != basicConst->get_1()) {
        auto temp_0 = add(n, basicConst->get_1());
        auto temp_1 = subtract(n, basicConst->get_1());
        *temp_1 /= *temp_0;
        auto copy_temp_1 = new Numerical(temp_1);
        auto temp_2 = getOne();

        copy_temp_1->a = temp_1->a = 0;
        while(true) {
            //Calculate one term of the taylor series.
            auto temp = *temp_1 / *temp_2;
            temp->a = 0;
            *result += *temp;
            delete temp;
            //Here the temp means the criteria of break.
            *temp_1 *= *copy_temp_1;
            *temp_2 += basicConst->get_1();
            temp = *temp_1 / *temp_2;
            int temp_power = temp->power;
            delete temp;

            if(result->power - temp_power >= basicConst->getGlobalPrecision())
                break;
            //Prepare for next calculate.
            *temp_1 *= *copy_temp_1;
            *temp_2 += basicConst->get_1();
        }
        *result *= basicConst->get_2();
        delete temp_0;
        delete temp_1;
        delete copy_temp_1;
        delete temp_2;
    }
    return result;
}

Numerical* ln(const Numerical& n) {
    auto result = ln_light(n);
    if(n.a != 0) {
        auto n_error = n.getMinimum();
        auto error = ln_light(*n_error);
        if(error == nullptr) {
            delete n_error;
            n_error = n.getMaximum();
            error = ln_light(*n_error);
        }
        *error -= *result;
        error->length = (signed char)error->getSize();

        result->applyError(error);

        delete n_error;
        delete error;
    }
    return result;
}
//Return log_a n
Numerical* log(const Numerical& n, const Numerical& a) {
    if(a == basicConst->get_1())
        return nullptr;

    auto ln_n = ln(n);
    if(ln_n == nullptr)
        return nullptr;

    auto ln_a = ln(a);
    if(ln_a == nullptr) {
        delete ln_n;
        return nullptr;
    }
    *ln_n /= *ln_a;
    delete ln_a;
    return ln_n;
}

Numerical* exp(const Numerical& n) {
    auto result = getOne();
    auto temp = new Numerical(n);
    auto rank = getOne();
    while(true) {
        *temp /= *rank;
        if(*temp < basicConst->getExpectedRelativeError())
            break;
        *result += *temp;
        *temp *= n;
        *rank += basicConst->get_1();
    }
    return result;
}

Numerical* pow(const Numerical& n, const Numerical& a) {
    auto result = ln(n);
    *result *= a;
    *result << *exp(*result);
    return result;
}
/*
 * Taylor's formula n.th term: (-1)^n * x^(2n) / (2n!)
 * Here temp_1 = x^(2n), temp_2 = 2n!, rank = 2n
 */
Numerical* cos(const Numerical& n) {
    auto result = getOne();
    if(n != basicConst->get_0()) {
        auto MachinePrecision = basicConst->getGlobalPrecision();
        auto& ONE = basicConst->get_1();

        auto square_n = n * n;
        auto temp_1 = new Numerical(square_n);
        auto temp_2 = getTwo();
        auto rank = getTwo();

        while(true) {
            //Calculate one term of the taylor series.
            auto temp = *temp_1 / *temp_2;
            temp->length = (signed char)-temp->length;
            *result += *temp;
            //Here the temp means the criteria of break.
            *temp *= n;
            *rank += ONE;
            *temp /= *rank;
            int temp_power = temp->power;
            delete temp;

            if(result->power - temp_power >= MachinePrecision)
                break;
            //Prepare for next calculate.
            *temp_1 *= *square_n;
            *temp_2 *= *rank;
            *rank += ONE;
            *temp_2 *= *rank;
        }
        delete square_n;
        delete temp_1;
        delete temp_2;
        delete rank;
    }
    return result;
}

Numerical* sin(const Numerical& n) {
    auto result = getZero();
    if(n != basicConst->get_0()) {
        auto MachinePrecision = basicConst->getGlobalPrecision();
        auto& ONE = basicConst->get_1();

        auto square_n = n * n;
        auto temp_1 = new Numerical(n);
        auto temp_2 = getOne();
        auto rank = getOne();

        while(true) {
            //Calculate one term of the taylor series.
            auto temp = *temp_1 / *temp_2;
            temp->length = (signed char)-temp->length;
            *result += *temp;
            //Here the temp means the criteria of break.
            *temp *= n;
            *rank += ONE;
            *temp /= *rank;
            int temp_power = temp->power;
            delete temp;

            if(result->power - temp_power >= MachinePrecision)
                break;
            //Prepare for next calculate.
            *temp_1 *= *square_n;
            *temp_2 *= *rank;
            *rank += ONE;
            *temp_2 *= *rank;
        }
        delete square_n;
        delete temp_1;
        delete temp_2;
        delete rank;
    }
    return result;
}

Numerical* tan(const Numerical& n) {
    auto cos_n = cos(n);
    if(cos_n == nullptr)
        return nullptr;
    auto sin_n = sin(n);
    *sin_n /= *cos_n;
    delete cos_n;
    return sin_n;
}

Numerical* sec(const Numerical& n) {
    auto cos_n = cos(n);
    auto result = reciprocal(*cos_n);
    delete cos_n;
    return result;
}

Numerical* csc(const Numerical& n) {
    auto sin_n = sin(n);
    auto result = reciprocal(*sin_n);
    delete sin_n;
    return result;
}

Numerical* cot(const Numerical& n) {
    auto sin_n = sin(n);
    if(sin_n == nullptr)
        return nullptr;
    auto cos_n = cos(n);
    *cos_n /= *sin_n;
    delete sin_n;
    return cos_n;
}

Numerical* arccos(const Numerical& n) {
    return bisectionMethod(cos, n, basicConst->get_0(), mathConst->getPI(), basicConst->get_1(), basicConst->getMinus_1());
}

Numerical* arcsin(const Numerical& n) {
    return bisectionMethod(sin, n, mathConst->getMinus_PI_2(), mathConst->getPI_2(), basicConst->getMinus_1(), basicConst->get_1());
}

Numerical* arctan(const Numerical& n) {
    auto temp = n * n;
    *temp += basicConst->get_1();
    auto sqrt_temp = sqrt(*temp);
    delete temp;

    temp = n / *sqrt_temp;
    delete sqrt_temp;

    auto result = arcsin(*temp);
    if((result->length ^ n.length) < 0) // NOLINT(hicpp-signed-bitwise)
        result->length = (signed char)-result->length;

    delete temp;
    return result;
}

Numerical* arcsec(const Numerical& n) {
    auto temp = reciprocal(n);
    if(temp == nullptr)
        return temp;
    auto result = arccos(*temp);
    delete temp;
    return result;
}

Numerical* arccsc(const Numerical& n) {
    auto temp = reciprocal(n);
    if(temp == nullptr)
        return temp;
    auto result = arcsin(*temp);
    delete temp;
    return result;
}

Numerical* arccot(const Numerical& n) {
    auto temp = reciprocal(n);
    if(temp == nullptr)
        return temp;
    auto result = arctan(*temp);
    delete temp;
    return result;
}

Numerical* cosh(const Numerical& n) {
    auto result = exp(n);
    auto temp = reciprocal(*result);
    *result += *temp;
    *result /= basicConst->get_2();
    delete temp;
    return result;
}

Numerical* sinh(const Numerical& n) {
    auto result = exp(n);
    auto temp = reciprocal(*result);
    *result -= *temp;
    *result /= basicConst->get_2();
    delete temp;
    return result;
}

Numerical* tanh(const Numerical& n) {
    auto result = exp(n);
    auto temp = reciprocal(*result);
    auto temp1 = *result + *temp;
    *result -= *temp;
    *result /= *temp1;
    delete temp;
    delete temp1;
    return result;
}

Numerical* sech(const Numerical& n) {
    auto result = getTwo();
    auto temp = exp(n);
    auto temp1 = reciprocal(*temp);
    *temp += *temp1;
    *result /= *temp;
    delete temp;
    delete temp1;
    return result;
}

Numerical* csch(const Numerical& n) {
    auto result = getTwo();
    auto temp = exp(n);
    auto temp1 = reciprocal(*temp);
    *temp -= *temp1;
    *result /= *temp;
    delete temp;
    delete temp1;
    return result;
}

Numerical* coth(const Numerical& n) {
    auto result = exp(n);
    auto temp = reciprocal(*result);
    auto temp1 = *result - *temp;
    *result += *temp;
    *result /= *temp1;
    delete temp;
    delete temp1;
    return result;
}

Numerical* arccosh(const Numerical& n) {
    auto temp = n * n;
    *temp -= basicConst->get_1();
    *temp << *sqrt(*temp);
    *temp += n;
    auto result = ln(*temp);
    delete temp;
    return result;
}

Numerical* arcsinh(const Numerical& n) {
    auto temp = n * n;
    *temp += basicConst->get_1();
    *temp << *sqrt(*temp);
    *temp += n;
    auto result = ln(*temp);
    delete temp;
    return result;
}

Numerical* arctanh(const Numerical& n) {
    auto result = basicConst->get_1() + n;
    auto temp = basicConst->get_1() - n;
    *result /= *temp;
    *result << *ln(*result);
    *result /= basicConst->get_2();
    delete temp;
    return result;
}

Numerical* arcsech(const Numerical& n) {
    auto temp = reciprocal(n);
    auto result = arccosh(*temp);
    delete temp;
    return result;
}

Numerical* arccsch(const Numerical& n) {
    auto temp = reciprocal(n);
    auto result = arcsinh(*temp);
    delete temp;
    return result;
}

Numerical* arccoth(const Numerical& n) {
    auto result = n + basicConst->get_1();
    auto temp = n - basicConst->get_1();
    *result /= *temp;
    *result << *ln(*result);
    *result /= basicConst->get_2();
    delete temp;
    return result;
}