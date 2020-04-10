/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 * This file contains some low levels which will be used by Numerical.
 */
#include <climits>
#include <cstring>
#include <QtCore/qlogging.h>
#include "CalcBasic.h"

#pragma clang diagnostic push
#pragma ide diagnostic ignored "hicpp-signed-bitwise"
extern const unsigned char ByteMask = -1; //1111 1111 in binary form.
extern const unsigned char ByteLowMask = ByteMask >> (CHAR_BIT / 2); //0000 1111 in binary form.
extern const unsigned char UnitByte = sizeof(char);

//n1 * n2 = product(16 bits) = carry(high 8 bits) + ReturnValue(low 8bits)
unsigned char basicMultiply(unsigned char& carry, unsigned char n1, unsigned char n2) {
    unsigned char n1_low = n1 & ByteLowMask;
    unsigned char n1_high = n1 >> (CHAR_BIT / 2);
    unsigned char n2_low = n2 & ByteLowMask;
    unsigned char n2_high = n2 >> (CHAR_BIT / 2);
    
    auto ll = n1_low * n2_low;
    auto lh = n1_low * n2_high;
    auto hl = n1_high * n2_low;
    auto hh = n1_high * n2_high;
    
    lh += ll >> (CHAR_BIT / 2);
    lh += hl;
    if(lh < hl)
        hh += (unsigned char)1 << (CHAR_BIT / 2);
    carry = hh + (lh >> (CHAR_BIT / 2));
    return (lh << (CHAR_BIT / 2)) + (ll & ByteLowMask);
}
#pragma clang diagnostic pop
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
            shallow_copy = new Numerical(n2.byte, (signed char)-n2.length, n2.power);
            result = n1 - *shallow_copy;
            shallow_copy->byte = nullptr;
        }
        else {
            shallow_copy = new Numerical(n1.byte, (signed char)-n1.length, n1.power);
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
        const int bigSize = big->getSize();
        const int smallSize = small->getSize();
        int lastIndex = smallSize - 1;
        //Estimate the ed of result first, will calculate it accurately later.
        signed char length = (signed char)(big->power + std::max(bigSize - big->power, smallSize - small->power));
        auto byte = (unsigned char*)malloc(length * sizeof(char));
        memcpy(byte + length - bigSize, big->byte, bigSize * sizeof(char));
        memset(byte, 0, (length - bigSize) * sizeof(char));

        int index = length - big->power + small->power - smallSize;
        unsigned char aByte;
        unsigned char carry = 0;
        unsigned char carry_temp;
        //Add small to big
        for(int i = 0; i < lastIndex; ++i) {
            aByte = byte[index];
            byte[index] += small->byte[i];
            carry_temp = aByte > byte[index];
            aByte = byte[index];
            byte[index] += carry;
            carry = carry_temp | (aByte > byte[index]);
            ++index;
        }
        aByte = byte[index];
        byte[index] += small->byte[lastIndex];
        carry_temp = aByte > byte[index];
        aByte = byte[index];
        byte[index] += carry;
        carry = carry_temp | (aByte > byte[index]);

        while(carry != 0 && index != length - 1) {
            ++index;
            aByte = byte[index];
            byte[index] += carry;
            carry = aByte > byte[index];
        }
        ///////////////////////////////////////Get byte, length and power//////////////////////////
        int power = big->power;
        if(carry != 0) {
            ++length;
            ++power;
            byte = (unsigned char*)realloc(byte, length * sizeof(char));
            byte[length - 1] = 1;
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
            shallow_copy = new Numerical(n2.byte, (signed char)-n2.length, n2.power);
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
redo:
            const int bigSize = big->getSize();
            const int smallSize = small->getSize();
            const int lastIndex = smallSize - 1;
            //Estimate the ed of result first, will calculate it accurately later.
            signed char length = (signed char)(big->power + std::max(bigSize - big->power, smallSize - small->power));
            auto byte = (unsigned char*)malloc(length * sizeof(char));
            memcpy(byte + length - bigSize, big->byte, bigSize * sizeof(char));
            memset(byte, 0, (length - bigSize) * sizeof(char));

            int index = length - big->power + small->power - smallSize;
            unsigned char aByte;
            unsigned char carry = 0;
            unsigned char carry_temp;
            //Subtract small from big
            for(int i = 0; i < lastIndex; ++i) {
                aByte = byte[index];
                byte[index] -= small->byte[i];
                carry_temp = aByte < byte[index];
                aByte = byte[index];
                byte[index] -= carry;
                carry = carry_temp | (aByte < byte[index]);
                ++index;
            }
            aByte = byte[index];
            byte[index] -= small->byte[lastIndex];
            carry_temp = aByte < byte[index];
            aByte = byte[index];
            byte[index] -= carry;
            carry = carry_temp | (aByte < byte[index]);

            while(carry != 0 && index != length - 1) {
                ++index;
                aByte = byte[index];
                byte[index] -= carry;
                carry = aByte < byte[index];
            }

            if(carry != 0) {
                auto temp = big;
                big = small;
                small = temp;
                changeSign = !changeSign;
                goto redo;
            }

            if(changeSign)
                length = (signed char)-length;
            result = new Numerical(byte, length, big->power);
            cutZero(result);
        }
    }
    else {
        shallow_copy = new Numerical(n1.byte, (signed char)-n1.length, n1.power);
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
        const int size1 = n1.getSize();
        const int size2 = n2.getSize();
        const int last2 = size2 - 1;
        //Estimate the ed of result first. we will calculate it accurately later.
        auto length = (signed char)(size1 + size2 - 1);
        auto byte = (unsigned char*)calloc(length, sizeof(char));
        unsigned char aByte;
        unsigned char carry = 0;
        //Every time the outer loop finished once, the position of carry will be reset. So we have to save the data.
        unsigned char carry_last = 0;
        unsigned char carry_temp;
        for (int i = 0; i < size1; ++i) {
            int index = i;
            for(int j = 0; j < last2; ++j) {
                aByte = byte[index];
                byte[index] += basicMultiply(carry_temp, n1.byte[i], n2.byte[j]);
                carry_temp += aByte > byte[index];
                aByte = byte[index];
                byte[index] += carry;
                carry = carry_temp + (aByte > byte[index]);
                ++index;
            }
            aByte = byte[index];
            byte[index] += basicMultiply(carry_temp, n1.byte[i], n2.byte[last2]);
            carry_temp += aByte > byte[index];
            aByte = byte[index];
            byte[index] += carry;
            carry_temp += aByte > byte[index];
            aByte = byte[index];
            byte[index] += carry_last;
            carry_last = carry_temp + (aByte > byte[index]);
        }
        ///////////////////////////////////////Get byte, length and power//////////////////////////;
        int power = n1.power + n2.power;
        if (carry_last != 0) {
            ++length;
            ++power;
            byte = (unsigned char*)realloc(byte, length * sizeof(char));
            byte[length - 1] = carry_last;
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
        //let n1_copy's power equal to n2_copy, power of the result will change correspondingly.
        int power = n1_copy->power - n2_copy->power;
        n1_copy->power = n2_copy->power;
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
        if((n1.length ^ n2.length) < 0) // NOLINT(hicpp-signed-bitwise)
            length = (signed char)-length;
        //1 comes from the algorithm
        result = new Numerical(byte, length, power, 1);
        cutZero(result);
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