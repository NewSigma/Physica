#include "../../Header/RealNumber.h"
#include <iostream>
#include <cstring>
#include <cmath>
#include "../../Header/Const.h"
#include "../../Header/BasicCalculates.h"

/*
 *
 * Useful formulas:
 * (1) n.ed - n.power
 *     Number of digits between the unit and the last digit. (Both included)
 *
 * Operators that do not need to free memory : = += -= *= /= > <
 * Operators that need memory free : new + - * / toRealNumber() randomRealNumber()
 *
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
extern const Const_1* const_1;
////////////////////////////////RealNumber////////////////////////////////
RealNumber::RealNumber() {
    byte = nullptr;
    length = power = a = 0;
    sign = true;
}
/*
 * May not be very accurate.
 */
RealNumber::RealNumber(double d, unsigned char acc) {
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
    } while(copy_d != 0 && length <= const_1->MachinePrecision);

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
RealNumber::operator double() const {
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

RealNumber::RealNumber(const RealNumber& n) {
    length = n.length;
    power = n.power;
    sign = n.sign;
    byte = (unsigned char*)malloc(length * sizeof(char));
    memcpy(byte, n.byte, length * sizeof(char));
    a = n.a;
}

RealNumber::RealNumber(unsigned char* b, int len, int pow, bool s, unsigned char acc) {
    byte = b;
    length = len;
    power = pow;
    sign = s;
    a = acc;
}

RealNumber::RealNumber(const RealNumber* n) {
    length = n->length;
    power = n->power;
    sign = n->sign;
    byte = (unsigned char*)malloc(length * sizeof(char));
    memcpy(byte, n->byte, length * sizeof(char));
    a = n->a;
}

RealNumber::~RealNumber() {
    free(byte);
}

std::ostream& operator<<(std::ostream& os, const RealNumber& n) {
    if(!n.sign)
        os << '-';
    if(n.power < 0) {
        os << "0.";
        for(int i = -n.power - 1; i > 0; --i)
            os << '0';
        for(int i = 0; i < n.length; ++i)
            os << (int)n.byte[i];
    }
    else {
        if(n.length > n.power) {
            for(int i = 0; i < n.power + 1; ++i)
                os << (int)n.byte[i];
            os << '.';
            for(int i = n.power + 1; i < n.length; ++i)
                os << (int)n.byte[i];
        }
        else {
            for(int i = 0; i < n.length; ++i)
                os << (int)n.byte[i];
            for(int i = n.power - n.length + 1; i > 0; --i)
                os << '0';
        }
    }
    os << "\tLength = " << n.length << "\tPower = " << n.power << "\tAccuracy = ";

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
void RealNumber::operator<<(RealNumber& n) {
    free(byte);
    byte = n.byte;
    length = n.length;
    power = n.power;
    sign = n.sign;
    a = n.a;
    n.byte = nullptr;
    delete &n;
}

RealNumber& RealNumber::operator= (const RealNumber& n) {
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

//Return accuracy in class RealNumber.
RealNumber* RealNumber::getAccuracy() const {
    auto b = (unsigned char*)malloc(sizeof(char));
    b[0] = a;
    return new RealNumber(b, 1, power - length + 1);
}
//Return this + accuracy
RealNumber* RealNumber::getMaximum() const {
    auto acc = getAccuracy();
    auto result = add(this, acc);
    delete acc;
    return result;
}
//Return this - accuracy
RealNumber* RealNumber::getMinimum() const {
    auto acc = getAccuracy();
    auto result = subtract(this, acc);
    delete acc;
    return result;
}
//Add error to this and adjust this->length as well as this-> byte.
bool RealNumber::applyError(const RealNumber* error) {
    int temp = power - length + 1 - error->power;
    if(temp <= 0) {
        if(temp < 0) {
            auto b = (unsigned char*)malloc(sizeof(char));
            b[0] = a;
            auto error_1 = new RealNumber(b, 1, power - length + 1);
            *error_1 << *add(error_1, error);
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
        std::cout << "[RealNumber] Warn: Accumulated too many errors." << std::endl;
        power += 1 - length;
        length = 1;
        return true;
    }

    byte = (unsigned char*)realloc(byte, length * sizeof(char));
    return false;
}

RealNumber* operator+ (const RealNumber& n1, const RealNumber& n2) {
    auto result = add(&n1, &n2);
    result->a = cutLength(result);
    if(n1.a != 0 || n2.a != 0) {
        RealNumber* error;
        if(n1.a == 0)
            error = n2.getAccuracy();
        else if(n2.a == 0)
            error = n1.getAccuracy();
        else {
            auto n1_a = n1.getAccuracy();
            auto n2_a = n2.getAccuracy();
            error = add(n1_a, n2_a);
            delete n1_a;
            delete n2_a;
        }
        result->applyError(error);
        delete error;
    }
    return result;
}

RealNumber* operator- (const RealNumber& n1, const RealNumber& n2) {
    auto result = subtract(&n1, &n2);
    result->a = cutLength(result);
    if(n1.a != 0 || n2.a != 0) {
        RealNumber* error;
        if(n1.a == 0)
            error = n2.getAccuracy();
        else if(n2.a == 0)
            error = n1.getAccuracy();
        else {
            auto n1_a = n1.getAccuracy();
            auto n2_a = n2.getAccuracy();
            error = add(n1_a, n2_a);
            delete n1_a;
            delete n2_a;
        }
        result->applyError(error);
        delete error;
    }
    return result;
}

RealNumber* operator* (const RealNumber& n1, const RealNumber& n2) {
    auto result = multiply(&n1, &n2);
    result->a = cutLength(result);
    if(n1.a != 0 || n2.a != 0) {
        RealNumber* error;
        if(n1.a == 0) {
            auto n2_a = n2.getAccuracy();
            error = multiply(&n1, n2_a);
            delete n2_a;
        }
        else if(n2.a == 0) {
            auto n1_a = n1.getAccuracy();
            error = multiply(&n2, n1_a);
            delete n1_a;
        }
        else {
            auto n1_a = n1.getAccuracy();
            auto n2_a = n2.getAccuracy();
            auto temp_1 = multiply(&n1, n2_a);
            auto temp_2 = multiply(&n2, n1_a);
            auto temp_3 = multiply(n1_a, n2_a);
            auto temp_4 = add(temp_1, temp_2);
            error = add(temp_3, temp_4);
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

RealNumber* operator/ (const RealNumber& n1, const RealNumber& n2) {
    auto result = divide(&n1, &n2);
    if(result == nullptr)
        return result;
    result->a += cutLength(result);
    if(n1.a != 0 || n2.a != 0) {
        RealNumber* error;
        if(n1.a == 0) {
            auto n2_a = n2.getAccuracy();
            auto temp_1 = multiply(&n1, n2_a);
            auto temp_2 = subtract(&n2, n2_a);
            auto temp_3 = multiply(&n2, temp_2);
            error = divide(temp_1, temp_3);
            delete n2_a;
            delete temp_1;
            delete temp_2;
            delete temp_3;
        }
        else if(n2.a == 0) {
            auto n1_a = n1.getAccuracy();
            error = divide(n1_a, &n2);
            delete n1_a;
        }
        else {
            auto n1_a = n1.getAccuracy();
            auto n2_a = n2.getAccuracy();
            auto temp_1 = multiply(&n1, n2_a);
            auto temp_2 = multiply(&n2, n1_a);
            auto temp_3 = add(temp_1, temp_2);
            auto temp_4 = subtract(&n2, n2_a);
            auto temp_5 = multiply(&n2, temp_4);
            error = divide(temp_3, temp_5);
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
    return result;
}
/*
 * Warning: n1 can not be temp object.
 */
void operator+= (RealNumber& n1, const RealNumber& n2) {
    RealNumber* p_result = n1 + n2;
    free(n1.byte);
    n1.byte = p_result->byte;
    n1.length = p_result->length;
    n1.power = p_result->power;
    n1.sign = p_result->sign;
    n1.a = p_result->a;
    p_result->byte = nullptr;
    delete p_result;
}

void operator-= (RealNumber& n1, const RealNumber& n2) {
    RealNumber* p_result = n1 - n2;
    free(n1.byte);
    n1.byte = p_result->byte;
    n1.length = p_result->length;
    n1.power = p_result->power;
    n1.sign = p_result->sign;
    n1.a = p_result->a;
    p_result->byte = nullptr;
    delete p_result;
}

void operator*= (RealNumber& n1, const RealNumber& n2) {
    RealNumber* p_result = n1 * n2;
    free(n1.byte);
    n1.byte = p_result->byte;
    n1.length = p_result->length;
    n1.power = p_result->power;
    n1.sign = p_result->sign;
    n1.a = p_result->a;
    p_result->byte = nullptr;
    delete p_result;
}
/*
 * n2 mustn't be zero.
 */
void operator/= (RealNumber& n1, const RealNumber& n2) {
    RealNumber* p_result = n1 / n2;
    free(n1.byte);
    n1.byte = p_result->byte;
    n1.length = p_result->length;
    n1.power = p_result->power;
    n1.sign = p_result->sign;
    n1.a = p_result->a;
    p_result->byte = nullptr;
    delete p_result;
}

bool operator> (const RealNumber& n1, const RealNumber& n2) {
    //Judge from sign.
    if(n1.isPositive()) {
        if(n2.isZero() || n2.isNegative())
            return true;
    }
    else if(n1.isZero())
        return n2.isNegative();
    else {
        if(n2.isPositive() || n2.isZero())
            return false;
    }
    //If we cannot get a result, judge from power
    bool result;
    if(n1.power > n2.power)
        result = true;
    else if(n1.power < n2.power)
        result = false;
    else {
        //The only method left.
        RealNumber* p_result = n1 - n2;
        result = p_result->sign;
        delete p_result;
    }
    return result;
}

bool operator< (const RealNumber& n1, const RealNumber& n2) {
    //Judge from sign.
    if(n1.isPositive()) {
        if(n2.isZero() || n2.isNegative())
            return false;
    }
    else if(n1.isZero())
        return n2.isPositive();
    else {
        if(n2.isPositive() || n2.isZero())
            return true;
    }
    //If we cannot get a result, judge from power
    bool result;
    if(n1.power > n2.power)
        result = false;
    else if(n1.power < n2.power)
        result = true;
    else {
        //The only method left.
        RealNumber* p_result = n1 - n2;
        result = p_result->isNegative();
        delete p_result;
    }
    return result;
}

bool operator>= (const RealNumber& n1, const RealNumber& n2) {
    return !(n1 < n2);
}

bool operator<= (const RealNumber& n1, const RealNumber& n2) {
    return !(n1 > n2);
}

bool operator== (const RealNumber& n1, const RealNumber& n2) {
    if(n1.power != n2.power)
        return false;
    if(n1.sign != n2.sign)
        return false;
    if(n1.a != n2.a)
        return false;
    const RealNumber* longer;
    const RealNumber* shorter;
    if(n1.length > n2.length) {
        longer = &n1;
        shorter = &n2;
    }
    else {
        longer = &n2;
        shorter = &n1;
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

bool operator!= (const RealNumber& n1, const RealNumber& n2) {
    return !(n1 == n2);
}
/*
 * Using Newton's method
 * May return nullptr.
 */
RealNumber* operator^ (const RealNumber& n1, const RealNumber& n2) {
    RealNumber* result = nullptr;
    if(n1.isZero()) {
        if(!n2.isZero())
            result = const_1->getZero();
    }
    else if(n1.isPositive()) {
        if(n2.isInteger()) {
            auto n2_copy = new RealNumber(n2);

            result = const_1->getOne();
            if(n2.isNegative()) {
                while(*n2_copy != *const_1->ZERO) {
                    *n2_copy -= *const_1->ONE;
                    *result /= n1;
                }
            }
            else {
                while(*n2_copy != *const_1->ZERO) {
                    *n2_copy -= *const_1->ONE;
                    *result *= n1;
                }
            }
            delete n2_copy;
        }
        else {
            //Do not have to use RealNumberA
            auto temp_result = const_1->getOne();
            auto temp_1 = ln(&n1);
            *temp_1 *= n2;
            *temp_1 += *const_1->ONE;
            bool go_on;
            do {
                result = ln(temp_result);
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

void operator^= (RealNumber& n1, const RealNumber& n2) {
    RealNumber* p_result = n1 ^ n2;
    n1 = *p_result;
    delete p_result;
}

RealNumber* operator- (const RealNumber& n) {
    auto result = new RealNumber(n);
    result->sign = !result->sign;
    return result;
}
////////////////////////////////Helper functions/////////////////////////////////////
//Return a real number between 0 and 1.
RealNumber* randomRealNumber() {
    srand(clock());
    srand(random());
    auto result = new RealNumber((double)random());
    *result /= *const_1->R_MAX;

    return result;
}
//Return a real number lowerBound and upperBound.
RealNumber* randomRealNumber(RealNumber* lowerBound, RealNumber* upperBound) {
    RealNumber* random = randomRealNumber();
    auto result = *lowerBound - *upperBound;
    *random *= *random;
    *result += *lowerBound;
    delete random;

    return result;
}
//////////////////////////////Process functions////////////////////////////////////////
/*
 * The following four functions simply calculate the result while operator functions will
 * consider the accuracy.
 */
RealNumber* add (const RealNumber* n1, const RealNumber* n2) {
    RealNumber* result;
    if(n1->isZero())
        result = new RealNumber(n2);
    else if(n2->isZero())
        result = new RealNumber(n1);
    else if (n1->sign != n2->sign) {
        RealNumber* shallow_copy;
        if (n1->sign) {
            shallow_copy = new RealNumber(n2->byte, n2->length, n2->power, true);
            result = *n1 - *shallow_copy;
            shallow_copy->byte = nullptr;
        }
        else {
            shallow_copy = new RealNumber(n1->byte, n1->length, n1->power, true);
            result = *n2 - *shallow_copy;
            shallow_copy->byte = nullptr;
        }
        delete shallow_copy;
    }
    else {
        const RealNumber* big;
        const RealNumber* small;
        if (n1->power > n2->power) {
            big = n1;
            small = n2;
        }
        else {
            big = n2;
            small = n1;
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
        ///////////////////////////////////////Get byte, length and power//////////////////////////
        unsigned char* byte;
        int power = big->power;
        if (temp[0] > 9) {
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
        result = new RealNumber(byte, length, power, big->sign);
    }
    return result;
}

RealNumber* subtract (const RealNumber* n1, const RealNumber* n2) {
    RealNumber* result;
    ///////////////////////////////////Deal with Binderies////////////////////////////////
    RealNumber* shallow_copy = nullptr;
    if(n1->isZero())
        result = -*n2;
    else if(n2->isZero())
        result = new RealNumber(n1);
    else if (n1->sign) {
        if (!n2->sign) {
            shallow_copy = new RealNumber(n2->byte, n2->length, n2->power, true);
            result = *n1 + *shallow_copy;
            shallow_copy->byte = nullptr;
        }
        else {
            ////////////////////////////////Calculate cursory first//////////////////////////////////////
            int unitIndex = std::max(n1->power, n2->power);
            //memcpy starts from the start index.
            int start = unitIndex - n1->power;
            //Estimate the ed of result first, will calculate it accurately later.
            int length = unitIndex + std::max(n1->length - n1->power, n2->length - n2->power);
            auto byte = (char*)malloc(length * sizeof(char));
            memset(byte, 0, start * sizeof(char));
            memcpy(byte + start, n1->byte, n1->length * sizeof(char));
            memset(byte + start + n1->length, 0, (length - start - n1->length) * sizeof(char));

            for (int i = n2->length - 1; i >= 0; --i) {
                int index = unitIndex - n2->power + i;
                byte[index] = char(byte[index] - n2->byte[i]);
            }

            for(int i = length - 1; i > 0; --i) {
                //Carry bit.
                if (byte[i] < 0) {
                    --byte[i - 1];
                    byte[i] += 10;
                }
            }
            //If n1 - n2 < 0, we have to change our method.
            if(byte[0] < 0) {
                free(byte);
                result = *n2 - *n1;
                result->sign = false;
            }
            else
                result = new RealNumber((unsigned char*)byte, length, unitIndex, true);
            cutZero(result);
        }
    }
    else {
        shallow_copy = new RealNumber(n1->byte, n1->length, n1->power, true);
        if (n2->sign) {
            result = *shallow_copy + *n2;
            result->sign = false;
        }
        else {
            auto shallow_copy_1 = new RealNumber(n2->byte, n2->length, n2->power, true);
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

RealNumber* multiply (const RealNumber* n1, const RealNumber* n2) {
    RealNumber* result;
    if(*n1 == *const_1->ONE)
        result = new RealNumber(n2);
    else if(*n2 == *const_1->ONE)
        result = new RealNumber(n1);
    else {
        //In this case, big has more digits after the dot.
        const RealNumber* big;
        const RealNumber* small;
        if (n1->length - n1->power > n2->length - n2->power) {
            big = n1;
            small = n2;
        }
        else {
            big = n2;
            small = n1;
        }
        ////////////////////////////////Calculate cursory first//////////////////////////////////////
        //Estimate the ed of result first. we will calculate it accurately later.
        int length = n1->length + n2->length - 1;
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
        int power = n1->power + n2->power;
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
        result = new RealNumber(byte, length, power, n1->sign == n2->sign);
    }
    return result;
}

RealNumber* divide (const RealNumber* n1, const RealNumber* n2) {
    if(n2->byte[0] == 0)
        return nullptr;

    RealNumber* result;
    if(n1->isZero() || n2->isZero())
        result = const_1->getZero();
    else if(*n2 == *const_1->ONE)
        result = new RealNumber(n1);
    else {
        result = new RealNumber();
        auto n1_copy = new RealNumber(n1);
        auto n2_copy = new RealNumber(n2);
        n1_copy->sign = n2_copy->sign = true;
        ////////////////////////////////Calculate cursory first//////////////////////////////////////
        //Estimate the ed of result first, we will calculate it accurately later.
        int length = const_1->MachinePrecision;
        int power = n1_copy->power - n2_copy->power - 1;
        auto temp = (unsigned char*)malloc(length * sizeof(char));
        memset(temp, 0, length * sizeof(char));

        auto n1_copy_old = n1_copy;
        n1_copy->power = n2_copy->power + 1;
        for (int i = 0; i < length; ++i) {
            char unit = 0;
            while(true) {
                n1_copy = subtract(n1_copy, n2_copy);
                if(n1_copy->isNegative()) {
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
        result->sign = n1->sign == n2->sign;
        cutZero(result);
    }
    return result;
}
/*
 * If the length of new array is larger than MachinePrecision, it will be set to MachinePrecision.
 * Return true array is cut.
 */
bool cutLength(RealNumber* n) {;
    bool result = false;
    int lastCutIndex = n->length;

    if(n->length > const_1->MachinePrecision) {
        lastCutIndex = const_1->MachinePrecision;
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
void cutZero(RealNumber* n) {
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

        if(n->byte[0] == 0)
            n->power = 0;
        else
            n->power -= firstCutIndex;
    }
}