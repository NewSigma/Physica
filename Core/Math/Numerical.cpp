/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "../Header/Numerical.h"
#include "../Header/Solve.h"
#include <cstring>
#include <cmath>

extern const Const_2* const_2;
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
Numerical::Numerical() {
    byte = nullptr;
    length = power = a = 0;
    sign = true;
}

Numerical::Numerical(std::wstring s, unsigned char acc) {
    byte = new unsigned char[s.size()];
    sign = true;
    a = acc;
    int index = 0, id_byte = 0, start_eff_num = 0, point_id = s.size();
    for(; index < s.size(); ++index) {
        switch(s[index]) {
            case '-':
                sign = false;
                break;
            case '.':
                point_id = index;
            case '0':
                break;
            default:
                start_eff_num = index;
                goto double_break;
        }
    }
    double_break:
    for(; index < s.size(); ++index) {
        switch(s[index]) {
            case '.':
                point_id = index;
                break;
            default:
                byte[id_byte] = s[index] - 48;
                ++id_byte;
        }
    }
    length = id_byte > const_1->GlobalPrecision ? const_1->GlobalPrecision : id_byte;
    power = point_id - start_eff_num + (point_id > start_eff_num ? -1 : 0);
    byte = (unsigned char*)realloc(byte, length);
}
/*
 * May not be very accurate.
 */
Numerical::Numerical(double d, unsigned char acc) {
    sign = d >= 0;
    d = sign ? d : -d;
    auto copy_d = d;

    length = power = 0;
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
    } while(copy_d != 0 && length <= const_1->GlobalPrecision);

    byte = (unsigned char*)malloc(length * sizeof(char));
    for(int i = 0; i < length; ++i) {
        byte[i] = (char)d;
        d -= double(int(d));
        d *= 10;
    }

    a = acc;
}
/*
 * Not very accurate either.
 */
Numerical::operator double() const {
    double result_integer = 0;
    double result_float = 0;

    int temp_index = power + 1;
    for(int i = 0; i < temp_index; ++i) {
        result_integer *= 10;
        result_integer += byte[i];
    }

    while(temp_index < length) {
        result_float += byte[temp_index];
        result_float *= 10;
        ++temp_index;
    }
    result_float /= pow(10,length- power);

    return result_integer + result_float;
}

Numerical::Numerical(const Numerical& n) {
    length = n.length;
    power = n.power;
    sign = n.sign;
    byte = (unsigned char*)malloc(length * sizeof(char));
    memcpy(byte, n.byte, length * sizeof(char));
    a = n.a;
}

Numerical::Numerical(unsigned char* b, int len, int pow, bool s, unsigned char acc) {
    byte = b;
    length = len;
    power = pow;
    sign = s;
    a = acc;
}

Numerical::Numerical(const Numerical* n) {
    length = n->length;
    power = n->power;
    sign = n->sign;
    byte = (unsigned char*)malloc(length * sizeof(char));
    memcpy(byte, n->byte, length * sizeof(char));
    a = n->a;
}

Numerical::~Numerical() {
    free(byte);
}

std::string Numerical::toString() const {
    std::string result;
    if(byte != nullptr) {
        if(isNegative())
            result.push_back('-');
        if(abs(power) > const_1->GlobalPrecision) {
            result.push_back(byte[0] + '0');
            result.push_back('.');
            for(int i = 1; i < length; ++i)
                result.push_back(byte[i] + '0');
            result += "×10^";
            result += std::to_string(power);
        }
        else if(power < 0) {
            result += "0.";
            for(int i = -power - 1; i > 0; --i)
                result.push_back('0');
            for(int i = 0; i < length; ++i)
                result.push_back(byte[i] + '0');
        }
        else {
            if(length > power) {
                for (int i = 0; i < power + 1; ++i)
                    result.push_back(byte[i] + '0');
                result.push_back('.');
                for(int i = power + 1; i < length; ++i)
                    result.push_back(byte[i] + '0');
            }
            else {
                for(int i = 0; i < length; ++i)
                    result.push_back(byte[i] + '0');
                for(int i = power - length + 1; i > 0; --i)
                    result.push_back('0');
            }
        }
    }
    return result;
}

std::ostream& operator<<(std::ostream& os, const Numerical& n) {
    os << n.toString() << "\tLength = " << n.length << "\tPower = " << n.power << "\tAccuracy = ";

    int temp = n.power - n.length +1;
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
    sign = n.sign;
    a = n.a;
    n.byte = nullptr;
    delete &n;
}

Numerical& Numerical::operator= (const Numerical& n) {
    if(this == &n)
        return *this;
    byte = (unsigned char*)realloc(byte, n.length * sizeof(char));
    memcpy(byte, n.byte, n.length * sizeof(char));
    length = n.length;
    power = n.power;
    sign = n.sign;
    a = n.a;
    return *this;
}

//Return accuracy in class Numerical.
Numerical* Numerical::getAccuracy() const {
    auto b = (unsigned char*)malloc(sizeof(char));
    b[0] = a;
    return new Numerical(b, 1, power - length + 1);
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
    int temp = power - length + 1 - error->power;
    if(temp <= 0) {
        if(temp < 0) {
            auto b = (unsigned char*)malloc(sizeof(char));
            b[0] = a;
            auto error_1 = new Numerical(b, 1, power - length + 1);
            *error_1 << *add(*error_1, *error);
            length += temp;
            a += error_1->byte[0];
            delete error_1;
        }
        else
            a += error->byte[0];
    }

    if(a > 9) {
        a = 1;
        --length;
    }

    if(length < 1) {
        std::cout << "[Numerical] Warn: Accumulated too many errors.\n";
        power += 1 - length;
        length = 1;
        return true;
    }

    byte = (unsigned char*)realloc(byte, length * sizeof(char));
    return false;
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
        result->a += cutLength(result);
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
    if(__glibc_unlikely(isZero())) {
        if(!n.isZero())
            result = getZero();
    }
    else if(isPositive()) {
        if(n.isInteger()) {
            auto n2_copy = new Numerical(n);

            result = getOne();
            if(n.isNegative()) {
                auto temp = reciprocal(*this);
                while(*n2_copy != *const_1->_0) {
                    *n2_copy -= *const_1->_1;
                    *result *= *temp;
                }
                delete temp;
            }
            else {
                while(*n2_copy != *const_1->_0) {
                    *n2_copy -= *const_1->_1;
                    *result *= *this;
                }
            }
            delete n2_copy;
        }
        else {
            auto temp_result = getOne();
            auto temp_1 = ln(*this);
            *temp_1 *= n;
            *temp_1 += *const_1->_1;
            bool go_on;
            do {
                result = ln(*temp_result);
                result->sign = !result->sign;
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
 * Warning: n1 can not be temp object.
 */
void Numerical::operator+= (const Numerical& n) {
    Numerical* p_result = *this + n;
    free(byte);
    byte = p_result->byte;
    length = p_result->length;
    power = p_result->power;
    sign = p_result->sign;
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
    sign = p_result->sign;
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
    sign = p_result->sign;
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
    sign = p_result->sign;
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
        result = p_result->sign;
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
    if(sign != n.sign)
        return false;
    if(a != n.a)
        return false;
    const Numerical* longer;
    const Numerical* shorter;
    if(length > n.length) {
        longer = this;
        shorter = &n;
    }
    else {
        longer = &n;
        shorter = this;
    }
    for(int i = 0; i < shorter->length; ++i) {
        if(shorter->byte[i] != longer->byte[i])
            return false;
    }
    //Extra parts should be zero.
    for(int i = shorter->length; i < longer->length; ++i) {
        if(longer->byte[i] != 0)
            return false;
    }
    return true;
}

bool Numerical::operator!= (const Numerical& n) const {
    return !(*this == n);
}

Numerical* Numerical::operator- () const {
    auto result = new Numerical(this);
    result->sign = !result->sign;
    return result;
}
////////////////////////////////Helper functions/////////////////////////////////////
Numerical* getZero() {
    return new Numerical(const_1->_0);
}

Numerical* getOne() {
    return new Numerical(const_1->_1);
}

Numerical* getTwo() {
    return new Numerical(const_1->_2);
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
    else if (n1.sign != n2.sign) {
        Numerical* shallow_copy;
        if (n1.sign) {
            shallow_copy = new Numerical(n2.byte, n2.length, n2.power, true);
            result = n1 - *shallow_copy;
            shallow_copy->byte = nullptr;
        }
        else {
            shallow_copy = new Numerical(n1.byte, n1.length, n1.power, true);
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
        ////////////////////////////////Calculate cursory first//////////////////////////////////////
        //Estimate the ed of result first, will calculate it accurately later.
        int length = big->power + std::max(big->length - big->power,small->length - small->power);
        auto temp = (unsigned char*)malloc(length * sizeof(char));
        memcpy(temp, big->byte, big->length * sizeof(char));
        memset(temp + big->length, 0, (length - big->length) * sizeof(char));
        for (int i = small->length - 1; i >= 0; --i) {
            int index = big->power - small->power + i;
            temp[index] += small->byte[i];
            //Carry bit.
            if (temp[index] > 9 && index > 0) {
                ++temp[index - 1];
                temp[index] -= 10;
            }
        }

        for(int i = big->power - small->power - 1; i > 0; --i) {
            if(temp[i] > 9) {
                ++temp[i - 1];
                temp[i] -= 10;
            }
        }
        ///////////////////////////////////////Get byte, length and power//////////////////////////
        unsigned char* byte;
        int power = big->power;
        if (__glibc_unlikely(temp[0] > 9)) {
            ++length;
            ++power;
            byte = (unsigned char*)malloc(length * sizeof(char));
            byte[0] = 1;
            byte[1] = temp[0] - 10;
            memcpy(byte + 2, temp + 1, (length - 2) * sizeof(char));
            free(temp);
        }
        else
            byte = (unsigned char*)realloc(temp, length * sizeof(char));
        ////////////////////////////////////Out put////////////////////////////////////////
        result = new Numerical(byte, length, power, big->sign);
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
    else if (n1.sign) {
        if (!n2.sign) {
            shallow_copy = new Numerical(n2.byte, n2.length, n2.power, true);
            result = n1 + *shallow_copy;
            shallow_copy->byte = nullptr;
        }
        else {
            ////////////////////////////////Calculate cursory first//////////////////////////////////////
            int unitIndex = std::max(n1.power, n2.power);
            //memcpy starts from the start index.
            int start = unitIndex - n1.power;
            //Estimate the ed of result first, will calculate it accurately later.
            int length = unitIndex + std::max(n1.length - n1.power, n2.length - n2.power);
            auto byte = (char*)malloc(length * sizeof(char));
            memset(byte, 0, start * sizeof(char));
            memcpy(byte + start, n1.byte, n1.length * sizeof(char));
            memset(byte + start + n1.length, 0, (length - start - n1.length) * sizeof(char));

            for (int i = n2.length - 1; i >= 0; --i) {
                int index = unitIndex - n2.power + i;
                byte[index] = char(byte[index] - n2.byte[i]);
            }

            for(int i = length - 1; i > 0; --i) {
                //Carry bit.
                if (byte[i] < 0) {
                    --byte[i - 1];
                    byte[i] += 10;
                }
            }
            //If n1 - n2 < 0, we have to change our method.
            if(__glibc_unlikely(byte[0] < 0)) {
                free(byte);
                result = n2 - n1;
                result->sign = false;
            }
            else
                result = new Numerical((unsigned char*)byte, length, unitIndex, true);
            cutZero(result);
        }
    }
    else {
        shallow_copy = new Numerical(n1.byte, n1.length, n1.power, true);
        if (n2.sign) {
            result = *shallow_copy + n2;
            result->sign = false;
        }
        else {
            auto shallow_copy_1 = new Numerical(n2.byte, n2.length, n2.power, true);
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
    if(n1 == *const_1->_1)
        result = new Numerical(n2);
    else if(n2 == *const_1->_1)
        result = new Numerical(n1);
    else {
        //In this case, big has more digits after the dot.
        const Numerical* big;
        const Numerical* small;
        if (n1.length - n1.power > n2.length - n2.power) {
            big = &n1;
            small = &n2;
        }
        else {
            big = &n2;
            small = &n1;
        }
        ////////////////////////////////Calculate cursory first//////////////////////////////////////
        //Estimate the ed of result first. we will calculate it accurately later.
        int length = n1.length + n2.length - 1;
        auto temp = (unsigned char*)malloc(length * sizeof(char));
        memset(temp, 0, length * sizeof(char));
        for (int i = small->length - 1; i >= 0; --i) {
            for(int j=big->length - 1; j >= 0; --j) {
                int index = i + j;
                temp[index] += small->byte[i] * big->byte[j];
                //Carry bit.
                if (temp[index] > 9 && index > 0) {
                    int tens = temp[index] / 10;
                    temp[index - 1] += tens;
                    temp[index] -= tens * 10;
                }
            }
        }
        ///////////////////////////////////////Get byte, length and power//////////////////////////
        unsigned char* byte;
        int power = n1.power + n2.power;
        if (temp[0] > 9) {
            ++length;
            ++power;
            int tens = temp[0] / 10;
            byte = (unsigned char*)malloc(length * sizeof(char));
            byte[0] = tens;
            byte[1] = temp[0] - tens * 10;
            memcpy(byte + 2, temp + 1, (length - 2) * sizeof(char));
            free(temp);
        }
        else
            byte = (unsigned char*)realloc(temp, length * sizeof(char));
        ////////////////////////////////////Out put////////////////////////////////////////
        result = new Numerical(byte, length, power, n1.sign == n2.sign);
    }
    return result;
}

Numerical* divide (const Numerical& n1, const Numerical& n2) {
    if(n2.byte[0] == 0)
        return nullptr;

    Numerical* result;
    if(__glibc_unlikely(n1.isZero() || n2.isZero()))
        result = getZero();
    else if(__glibc_unlikely(n2 == *const_1->_1))
        result = new Numerical(n1);
    else {
        result = new Numerical();
        auto n1_copy = new Numerical(n1);
        auto n2_copy = new Numerical(n2);
        n1_copy->sign = n2_copy->sign = true;
        ////////////////////////////////Calculate cursory first//////////////////////////////////////
        //Estimate the ed of result first, we will calculate it accurately later.
        int length = const_1->GlobalPrecision;
        int power = n1_copy->power - n2_copy->power - 1;
        auto temp = (unsigned char*)malloc(length * sizeof(char));
        memset(temp, 0, length * sizeof(char));

        auto n1_copy_old = n1_copy;
        n1_copy->power = n2_copy->power + 1;
        for (int i = 0; i < length; ++i) {
            char unit = 0;
            while(true) {
                n1_copy = subtract(*n1_copy, *n2_copy);
                if(__glibc_unlikely(n1_copy->isNegative())) {
                    delete n1_copy;
                    n1_copy = n1_copy_old;
                    break;
                }
                else {
                    ++unit;
                    delete n1_copy_old;
                    n1_copy_old = n1_copy;
                    if(n1_copy->isZero()) {
                        temp[i] = unit;
                        //Do we have better choice?
                        goto double_break;
                    }
                }
            }
            ++n1_copy->power;
            temp[i] = unit;
        }
        //1 comes from algorithm of divide()
        result->a = 1;
        double_break:
        delete n1_copy;
        delete n2_copy;
        ////////////////////////////////////Out put////////////////////////////////////////
        unsigned char* byte;
        if(temp[0] > 9) {
            ++power;
            byte = (unsigned char*)malloc(length * sizeof(char));
            byte[0] = temp[0] / 10;
            byte[1] = temp[0] - byte[0] * 10;
            memcpy(byte + 2, temp + 1, (length - 2) * sizeof(char));
            free(temp);
        }
        else
            byte = temp;
        result->byte = byte;
        result->length = length;
        result->power = power;
        result->sign = n1.sign == n2.sign;
        cutZero(result);
    }
    return result;
}
/*
 * If the length of new array is larger than GlobalPrecision, it will be set to GlobalPrecision.
 * Return true array is cut.
 */
bool cutLength(Numerical* n) {
    bool result = false;
    int lastCutIndex = n->length;

    if(n->length > const_1->GlobalPrecision) {
        lastCutIndex = const_1->GlobalPrecision;
        result = true;
    }

    if(lastCutIndex != n->length) {
        n->length = lastCutIndex;
        n->byte = (unsigned char*)realloc(n->byte, lastCutIndex * sizeof(char));
    }
    return result;
}
/*
 * Cut zeros from the beginning.
 */
void cutZero(Numerical* n) {
    int firstCutIndex = 0;
    //Ignore zeros from the first index.
    while(n->byte[firstCutIndex] == 0 && firstCutIndex < n->length - 1)
        ++firstCutIndex;

    if(firstCutIndex != 0) {
        n->length -= firstCutIndex;

        auto new_array = (unsigned char*)malloc(n->length * sizeof(char));
        memcpy(new_array, n->byte + firstCutIndex, n->length * sizeof(char));
        free(n->byte);
        n->byte = new_array;

        if(__glibc_unlikely(n->byte[0] == 0))
            n->power = 0;
        else
            n->power -= firstCutIndex;
    }
}
//////////////////////////////Basic Operations////////////////////////////////////////
//Return a real number between 0 and 1.
Numerical* randomNumerical() {
    srand(clock());
    srand(random());
    auto result = new Numerical((double)random());
    *result /= *const_1->R_MAX;

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
    return *const_1->_1 / n;
}
/*
 * *_light functions do not consider the error caused by a. For example, sqrt_light does not calculate
 * (sqrt(n + a) - sqrt(n)) for error.
 */
Numerical* sqrt_light(const Numerical& n) {
    if(!n.sign)
        return nullptr;
    auto MachinePrecision = const_1->GlobalPrecision;
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
        *result /= *const_1->_2;
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
        error->sign = true;

        result->applyError(error);

        delete n_error;
        delete error;
    }
    return result;
}
//TODO not completed: Use gamma function.
Numerical* factorial(const Numerical& n) {
    if(!n.sign) {
        std::cout << "[BasicCalculates] Error: Cannot solve the factorial of a minus value." << std::endl;
        return nullptr;
    }

    Numerical* result;
    if(n.isInteger()) {
        result = getOne();
        auto temp = getOne();
        while(*temp < n) {
            *temp += *const_1->_1;
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
    if(n != *const_1->_1) {
        auto temp_0 = add(n, *const_1->_1);
        auto temp_1 = subtract(n, *const_1->_1);
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
            *temp_2 += *const_1->_1;
            temp = *temp_1 / *temp_2;
            int temp_power = temp->power;
            delete temp;

            if(result->power - temp_power >= const_1->GlobalPrecision)
                break;
            //Prepare for next calculate.
            *temp_1 *= *copy_temp_1;
            *temp_2 += *const_1->_1;
        }
        *result *= *const_1->_2;
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
        error->sign = true;

        result->applyError(error);

        delete n_error;
        delete error;
    }
    return result;
}
//Return log_a n
Numerical* log(const Numerical& n, const Numerical& a) {
    if(a == *const_1->_1)
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
        if(*temp < *const_1->expectedRelativeError)
            break;
        *result += *temp;
        *temp *= n;
        *rank += *const_1->_1;
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
    if(n != *const_1->_0) {
        auto MachinePrecision = const_1->GlobalPrecision;
        auto ONE = const_1->_1;

        auto square_n = n * n;
        auto temp_1 = new Numerical(square_n);
        auto temp_2 = getTwo();
        auto rank = getTwo();
        bool sign = false;

        while(true) {
            //Calculate one term of the taylor series.
            auto temp = *temp_1 / *temp_2;
            temp->sign = sign;
            *result += *temp;
            //Here the temp means the criteria of break.
            *temp *= n;
            *rank += *ONE;
            *temp /= *rank;
            int temp_power = temp->power;
            delete temp;

            if(result->power - temp_power >= MachinePrecision)
                break;
            //Prepare for next calculate.
            sign = !sign;
            *temp_1 *= *square_n;
            *temp_2 *= *rank;
            *rank += *ONE;
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
    if(n != *const_1->_0) {
        auto MachinePrecision = const_1->GlobalPrecision;
        auto ONE = const_1->_1;

        auto square_n = n * n;
        auto temp_1 = new Numerical(n);
        auto temp_2 = getOne();
        auto rank = getOne();
        bool sign = true;

        while(true) {
            //Calculate one term of the taylor series.
            auto temp = *temp_1 / *temp_2;
            temp->sign = sign;
            *result += *temp;
            //Here the temp means the criteria of break.
            *temp *= n;
            *rank += *ONE;
            *temp /= *rank;
            int temp_power = temp->power;
            delete temp;

            if(result->power - temp_power >= MachinePrecision)
                break;
            //Prepare for next calculate.
            sign = !sign;
            *temp_1 *= *square_n;
            *temp_2 *= *rank;
            *rank += *ONE;
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
    return bisectionMethod(cos, n, *const_1->_0, *const_2->PI, *const_1->_1, *const_1->Minus_1);
}

Numerical* arcsin(const Numerical& n) {
    return bisectionMethod(sin, n, *const_2->Minus_PI_2, *const_2->PI_2, *const_1->Minus_1, *const_1->_1);
}

Numerical* arctan(const Numerical& n) {
    auto temp = n * n;
    *temp += *const_1->_1;
    auto sqrt_temp = sqrt(*temp);
    delete temp;

    temp = n / *sqrt_temp;
    delete sqrt_temp;

    auto result = arcsin(*temp);
    result->sign = n.sign;
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
    *result /= *const_1->_2;
    delete temp;
    return result;
}

Numerical* sinh(const Numerical& n) {
    auto result = exp(n);
    auto temp = reciprocal(*result);
    *result -= *temp;
    *result /= *const_1->_2;
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
    *temp -= *const_1->_1;
    *temp << *sqrt(*temp);
    *temp += n;
    auto result = ln(*temp);
    delete temp;
    return result;
}

Numerical* arcsinh(const Numerical& n) {
    auto temp = n * n;
    *temp += *const_1->_1;
    *temp << *sqrt(*temp);
    *temp += n;
    auto result = ln(*temp);
    delete temp;
    return result;
}

Numerical* arctanh(const Numerical& n) {
    auto result = *const_1->_1 + n;
    auto temp = *const_1->_1 - n;
    *result /= *temp;
    *result << *ln(*result);
    *result /= *const_1->_2;
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
    auto result = n + *const_1->_1;
    auto temp = n - *const_1->_1;
    *result /= *temp;
    *result << *ln(*result);
    *result /= *const_1->_2;
    delete temp;
    return result;
}