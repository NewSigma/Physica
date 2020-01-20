#include "../../Header/RealNumber.h"
#include <iostream>
#include <cstring>
#include <cmath>
#include "../../Header/Const.h"
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
    length = power = 0;
    sign = true;
}
/*
 * May not be very accurate.
 */
RealNumber::RealNumber(double d) {
    sign = d >= 0;
    d = sign ? d : -d;

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

    auto copy_f = d;
    do {
        copy_f -= (float)((int)copy_f);
        copy_f *= 10;
        ++length;
    } while(copy_f != 0);
    length = length > const_1->MachinePrecision ? const_1->MachinePrecision : length;

    byte = new unsigned char[length];
    for(int i = 0; i < length; ++i) {
        byte[i] = (char)d;
        d -= (float)((int)d);
        d *= 10;
    }
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
    byte = new unsigned char[length];
    memcpy(byte, n.byte, length * sizeof(char));
}

RealNumber::RealNumber(unsigned char* b, int len, int pow, bool s) {
    byte = b;
    length = len;
    power = pow;
    sign = s;
}

RealNumber::RealNumber(const RealNumber* n) {
    length = n->length;
    power = n->power;
    sign = n->sign;
    byte = new unsigned char[length];
    memcpy(byte, n->byte, length * sizeof(char));
}

RealNumber::~RealNumber() {
    delete[] byte;
}

void RealNumber::print() const {
    std::cout << (double)(*this) << "\tLength = " << length << "\tPower = " << power;
}

RealNumber& RealNumber::operator= (const RealNumber& n) {
    if(this == &n)
        return *this;
    delete[] byte;
    byte = new unsigned char[n.length];
    memcpy(byte, n.byte, n.length * sizeof(char));
    length = n.length;
    power = n.power;
    sign = n.sign;
    return *this;
}

RealNumber* operator+ (const RealNumber& n1, const RealNumber& n2) {
    auto result = add(&n1, &n2);
    cutArray(result);
    return result;
}

RealNumber* operator- (const RealNumber& n1, const RealNumber& n2) {
    RealNumber* result = subtract(&n1, &n2);
    cutArray(result);
    return result;
}

RealNumber* operator* (const RealNumber& n1, const RealNumber& n2) {
    auto result = multiply(&n1, &n2);
    cutArray(result);
    return result;
}

RealNumber* operator/ (const RealNumber& n1, const RealNumber& n2) {
    auto result = divide(&n1, &n2);
    cutArray(result);
    return result;
}
/*
 * Warning: n1 can not be temp object.
 */
void operator+= (RealNumber& n1, const RealNumber& n2) {
    RealNumber* p_result = n1 + n2;
    n1 = *p_result;
    delete p_result;
}

void operator-= (RealNumber& n1, const RealNumber& n2) {
    RealNumber* p_result = n1 - n2;
    n1 = *p_result;
    delete p_result;
}

void operator*= (RealNumber& n1, const RealNumber& n2) {
    RealNumber* p_result = n1 * n2;
    n1 = *p_result;
    delete p_result;
}

void operator/= (RealNumber& n1, const RealNumber& n2) {
    RealNumber* p_result = n1 / n2;
    n1 = *p_result;
    delete p_result;
}

bool operator> (const RealNumber& n1, const RealNumber& n2) {
    RealNumber* p_result = n1 - n2;
    bool result = p_result->sign;
    delete p_result;
    return result;
}

bool operator< (const RealNumber& n1, const RealNumber& n2) {
    RealNumber* p_result = n1 - n2;
    bool result = !p_result->sign;
    delete p_result;
    return result;
}

bool operator>= (const RealNumber& n1, const RealNumber& n2) {
    RealNumber* p_result = n1 - n2;
    bool result = p_result->sign || p_result->byte[0] == 0;
    delete p_result;
    return result;
}

bool operator<= (const RealNumber& n1, const RealNumber& n2) {
    RealNumber* p_result = n1 - n2;
    bool result = !p_result->sign || p_result->byte[0] == 0;
    delete p_result;
    return result;
}

bool operator== (const RealNumber& n1, const RealNumber& n2) {
    RealNumber* p_result = n1 - n2;
    bool result = p_result->byte[0] == 0;
    delete p_result;
    return result;
}

bool operator!= (const RealNumber& n1, const RealNumber& n2) {
    RealNumber* p_result = n1 - n2;
    bool result = p_result->byte[0] != 0;
    delete p_result;
    return result;
}
//TODO not completed. n2 should be a float number
RealNumber* operator^ (const RealNumber& n1, const RealNumber& n2) {
    RealNumber* result = nullptr;
    if(isInteger(&n2)) {
        auto n1_copy = new RealNumber(n1);

        result = getOne();
        while(*n1_copy != *const_1->ZERO) {
            *n1_copy -= *const_1->ONE;
            *result *= n1;
        }
        delete n1_copy;
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
///////////////////////////////////////RealNumberA/////////////////////////////////////////
RealNumberA::RealNumberA() {
    byte = nullptr;
    length = power = 0;
    this->sign = true;
}

RealNumberA::RealNumberA(double d, char acc) : RealNumber(d) {
    a = acc;
}

RealNumberA::RealNumberA(const RealNumberA& n)  : RealNumber(n) {
    a = n.a;
}

RealNumberA::RealNumberA(unsigned char* byte, int length, int power, bool sign, unsigned char acc) : RealNumber(byte, length, power, sign) {
    a = acc;
}

RealNumberA::RealNumberA(const RealNumber* n, unsigned char acc) : RealNumber(n) {
    a = acc;
}

void RealNumberA::print() const {
    RealNumber::print();
    std::cout << "\tAccuracy = ";
    if(a == 0)
        std::cout << "0.";
    else {
        if(power > 0) {
            std::cout << (int)a;
            for(int i = 0; i < power; ++i)
                std::cout << "0";
            std::cout << ".";
        }
        else if(power == 0)
            std::cout << (int)a << ".";
        else {
            std::cout << "0.";
            for (int i = 0; i > power + 1; --i)
            {
                std::cout << "0";
            }
            std::cout << (int)a;
        }
    }
}

RealNumberA& RealNumberA::operator= (const RealNumberA& n) {
    if(this == &n)
        return *this;
    *(RealNumber*)this = (RealNumber)n;
    a = n.a;
    return *this;
}

RealNumberA* operator+ (const RealNumberA& n1, const RealNumberA& n2) {
    auto raw_result = add(&n1, &n2);
    unsigned char a = cutArray(raw_result);
    RealNumberA* result;
    if(n1.a == 0 && n2.a == 0)
        result = new RealNumberA(raw_result, a);
    else {
        //Get a
        int length = raw_result->length;
        if(n1.length - n1.power == n2.length - n2.power)
            a += n1.a + n2.a;
        else
            a += n1.length - n1.power > n2.length - n2.power ? n2.a + 1 : n1.a + 1;

        if(a > 9) {
            if(length == 1) {
                std::cout << "[RealNumber] Warn: Accumulated too many errors.";
                exit(EXIT_FAILURE);
            }
            a = 2;
            --length;
            if(raw_result->byte[raw_result->length - 1] >= 5)
                raw_result->byte[raw_result->length - 2] += 1;
        }
        //Get byte
        auto byte = new unsigned char[length];
        memcpy(byte, raw_result->byte, length * sizeof(char));
        int power = raw_result->power;
        bool sign = raw_result->sign;
        result = new RealNumberA(byte, length, power, sign, a);
    }
    delete raw_result;
    return result;
}

RealNumberA* operator- (const RealNumberA& n1, const RealNumberA& n2) {
    auto raw_result = subtract(&n1, &n2);
    unsigned char a = cutArray(raw_result);
    RealNumberA* result;
    if(n1.a == 0 && n2.a == 0)
        result = new RealNumberA(raw_result, a);
    else {
        //Get a
        int length = raw_result->length;
        if(n1.length - n1.power == n2.length - n2.power)
            a += n1.a + n2.a;
        else
            a += n1.length - n1.power > n2.length - n2.power ? (n2.a + 1) : (n1.a + 1);
        if(a > 9) {
            if(length == 1) {
                std::cout << "[RealNumber] Warn: Accumulated too many errors.";
                exit(EXIT_FAILURE);
            }
            a += 5;
            a = a / 10 + a % 10 != 0;
            --length;
            if(raw_result->byte[raw_result->length - 1] >= 5)
                raw_result->byte[raw_result->length - 2] += 1;
        }
        //Get byte
        auto byte = new unsigned char[length];
        memcpy(byte, raw_result->byte, length * sizeof(char));
        int power = raw_result->power;
        bool sign = raw_result->sign;
        result = new RealNumberA(byte, length, power, sign, a);
    }
    delete raw_result;
    return result;
}

RealNumberA* operator* (const RealNumberA& n1, const RealNumberA& n2) {
    auto raw_result = multiply(&n1, &n2);
    unsigned char a = cutArray(raw_result);
    RealNumberA* result;
    if(n1.a == 0 && n2.a == 0)
        result = new RealNumberA(raw_result, a);
    else {
        //Get a
        auto copy_byte_1 = new unsigned char[n1.length];
        memcpy(copy_byte_1, n1.byte, n1.length);
        auto copy_byte_2 = new unsigned char[n2.length];
        memcpy(copy_byte_2, n2.byte, n2.length);
        auto n1_a = RealNumber(copy_byte_1, 1, -(n1.length - n1.power - 1), true);
        auto n2_a = RealNumber(copy_byte_2, 1, -(n2.length - n2.power - 1), true);
        auto error = n1 * n2_a;
        auto error_1 = n2 * n1_a;
        auto error_2 = n1_a * n2_a;
        *error += *error_1;
        *error += *error_2;

        a += error->byte[0] + 1;
        int length = raw_result->length;
        if(a > 9) {
            if(length == 1) {
                std::cout << "[RealNumber] Warn: Accumulated too many errors.";
                exit(EXIT_FAILURE);
            }
            a = 2;
            --length;
        }
        delete error;
        delete error_1;
        delete error_2;
        //Get byte
        auto byte = new unsigned char[length];
        memcpy(byte, raw_result->byte, length * sizeof(char));
        int power = raw_result->power;
        bool sign = raw_result->sign;
        result = new RealNumberA(byte, length, power, sign, a);
    }
    delete raw_result;
    return result;
}

RealNumberA* operator/ (const RealNumberA& n1, const RealNumberA& n2) {
    auto raw_result = divide(&n1, &n2);
    unsigned char a = cutArray(raw_result);
    RealNumberA* result;
    if(n1.a == 0 && n2.a == 0)
        result = new RealNumberA(raw_result, a);
    else {
        //Get a
        RealNumber* error;
        auto n1_a = RealNumber(new unsigned char[1]{n1.a}, 1, -(n1.length - n1.power - 1), true);
        auto n2_a = RealNumber(new unsigned char[1]{n2.a}, 1, -(n2.length - n2.power - 1), true);
        auto numerator = n1 * n2_a;
        auto numerator_1 = n2 * n1_a;
        *numerator += *numerator_1;
        auto denominator = n2 - n2_a;
        *denominator *= n2;
        error = *numerator / *denominator;

        a += error->byte[0] + 1;
        int length = raw_result->length;
        if(a > 9) {
            if(length == 1) {
                std::cout << "[RealNumber] Warn: Accumulated too many errors.";
                exit(EXIT_FAILURE);
            }
            a = 2;
            --length;
        }
        delete numerator;
        delete numerator_1;
        delete denominator;
        delete error;
        //Get byte
        auto byte = new unsigned char[length];
        memcpy(byte, raw_result->byte, length * sizeof(char));
        int power = raw_result->power;
        bool sign = raw_result->sign;
        result = new RealNumberA(byte, length, power, sign, a);
    }
    delete raw_result;
    return result;
}

void operator+= (RealNumberA& n1, const RealNumberA& n2) {
    RealNumberA* p_result = n1 + n2;
    n1 = *p_result;
    delete p_result;
}

void operator-= (RealNumberA& n1, const RealNumberA& n2) {
    RealNumberA* p_result = n1 - n2;
    n1 = *p_result;
    delete p_result;
}

void operator*= (RealNumberA& n1, const RealNumberA& n2) {
    RealNumberA* p_result = n1 * n2;
    n1 = *p_result;
    delete p_result;
}

void operator/= (RealNumberA& n1, const RealNumberA& n2) {
    RealNumberA* p_result = n1 / n2;
    n1 = *p_result;
    delete p_result;
}

RealNumberA* operator- (const RealNumberA& n) {
    auto result = new RealNumberA(n);
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

bool isInteger(const RealNumber* n) {
    return n->length == n->power + 1;
}
//////////////////////////////Process functions////////////////////////////////////////
/*
 * The following four functions are public parts owned by RealNumber and RealNumberA.
 */
RealNumber* add (const RealNumber* n1, const RealNumber* n2) {
    RealNumber* result;
    if (n1->sign != n2->sign) {
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
        auto temp = new unsigned char[length];
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
            byte = new unsigned char[length];
            byte[0] = 1;
            byte[1] = temp[0] - 10;
            memcpy(byte + 2, temp + 1, (length - 2) * sizeof(char));
        }
        else {
            byte = new unsigned char[length];
            memcpy(byte, temp, length * sizeof(char));
        }
        delete[] temp;
        ////////////////////////////////////Out put////////////////////////////////////////
        result = new RealNumber(byte, length, power, big->sign);
    }
    return result;
}

RealNumber* subtract (const RealNumber* n1, const RealNumber* n2) {
    RealNumber* result;
    ///////////////////////////////////Deal with Binderies////////////////////////////////
    RealNumber* shallow_copy;
    if (n1->sign) {
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
            auto byte = new char[length];
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
                delete[] byte;
                result = *n2 - *n1;
                result->sign = false;
            }
            else
                result = new RealNumber((unsigned char*)byte, length, unitIndex, true);
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
        delete shallow_copy;
    }
    ////////////////////////////////////Out put////////////////////////////////////////
    return result;
}

RealNumber* multiply (const RealNumber* n1, const RealNumber* n2) {
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
    auto temp = new unsigned char[length]{0};
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
        byte = new unsigned char[length];
        byte[0] = tens;
        byte[1] = temp[0] - tens * 10;
        memcpy(byte + 2, temp + 1, (length - 2) * sizeof(char));
    }
    else {
        byte = new unsigned char[length];
        memcpy(byte, temp, length * sizeof(char));
    }
    delete[] temp;
    ////////////////////////////////////Out put////////////////////////////////////////
    return new RealNumber(byte, length, power, n1->sign == n2->sign);
}

RealNumber* divide (const RealNumber* n1, const RealNumber* n2) {
    if(n2->byte[0] == 0)
        throw std::invalid_argument("[RealNumber] Can not divide by zero!");

    auto n1_copy = new RealNumber(n1);
    auto n2_copy = new RealNumber(n2);
    n1_copy->sign = true;
    n2_copy->sign = true;
    ////////////////////////////////Calculate cursory first//////////////////////////////////////
    //Estimate the ed of result first, we will calculate it accurately later.
    /*
     * Here we add 1 to length, making it convenient to cutArray to judge whether the calculate is stopped
     * by the limitation of precision.
     */
    int length = const_1->MachinePrecision + 1;
    int power = n1_copy->power - n2_copy->power - 1;
    auto temp = new unsigned char[length];
    n1_copy->power = n2_copy->power + 1;
    for (int i = 0; i < length; ++i) {
        char unit = 0;
        while(true) {
            *n1_copy -= *n2_copy;
            if(n1_copy->sign)
                unit += 1;
            else {
                *n1_copy += *n2_copy;
                break;
            }
        }
        ++n1_copy->power;
        temp[i] = unit;
    }
    delete n1_copy;
    delete n2_copy;
    ////////////////////////////////////Out put////////////////////////////////////////
    unsigned char* byte;
    if(temp[0] > 9) {
        ++power;
        byte = new unsigned char[length];
        byte[0] = temp[0] / 10;
        byte[1] = temp[0] - byte[0] * 10;
        memcpy(byte + 2, temp + 1, (length - 2) * sizeof(char));
        delete[] temp;
    }
    else
        byte = temp;
    return new RealNumber(byte, length, power, n1->sign == n2->sign);
}
/*
 * Zero on both sides will be cut,
 * Elements whose indexes satisfy firstCutIndex <= index < lastCutIndex is not zero.
 * If the length of new array is larger than MachinePrecision, it will be set to MachinePrecision.
 *
 * Return true if the length of new array is larger than MachinePrecision and we will cut it.
 */
bool cutArray(RealNumber* n) {
    int MachinePrecision = const_1->MachinePrecision;
    bool result = false;
    unsigned char* byte = n->byte;
    int firstCutIndex = 0;
    int lastCutIndex = n->length;
    //Ignore zeros from the first index.
    while(byte[firstCutIndex] == 0 && firstCutIndex < lastCutIndex - 1)
        firstCutIndex += 1;

    int maxIndex = firstCutIndex + MachinePrecision;
    if(n->length > maxIndex) {
        lastCutIndex = maxIndex;
        result = true;
    }

    if(firstCutIndex != 0 || lastCutIndex != n->length) {
        n->length = lastCutIndex - firstCutIndex;

        auto new_array = new unsigned char[n->length];
        memcpy(new_array, byte + firstCutIndex, n->length * sizeof(char));
        delete[] byte;
        n->byte = new_array;

        if(n->byte[0] == 0)
            n->power = 0;
        else
            n->power -= firstCutIndex;
    }
    return result;
}