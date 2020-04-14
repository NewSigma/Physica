/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 * This file contains some low levels which will be used by Numerical.
 */
#include <cstring>
#include <QtCore/qlogging.h>
#include <climits>
#include "Numerical.h"

const unsigned long LongLowMask = ULONG_MAX >> (LONG_WIDTH / 2); // NOLINT(hicpp-signed-bitwise)

#pragma clang diagnostic push
#pragma ide diagnostic ignored "hicpp-signed-bitwise"
//n1 * n2 = product(16 bits) = carry(high 8 bits) + ReturnValue(low 8bits)
unsigned long basicMultiply(unsigned long& carry, unsigned long n1, unsigned long n2) {
    unsigned long n1_low = n1 & LongLowMask;
    unsigned long n1_high = n1 >> (LONG_WIDTH / 2);
    unsigned long n2_low = n2 & LongLowMask;
    unsigned long n2_high = n2 >> (LONG_WIDTH / 2);

    auto ll = n1_low * n2_low;
    auto lh = n1_low * n2_high;
    auto hl = n1_high * n2_low;
    auto hh = n1_high * n2_high;

    lh += ll >> (LONG_WIDTH / 2);
    lh += hl;
    if(lh < hl)
        hh += (unsigned long)1 << (LONG_WIDTH / 2);
    carry = hh + (lh >> (LONG_WIDTH / 2));
    return (lh << (LONG_WIDTH / 2)) + (ll & LongLowMask);
}
#pragma clang diagnostic pop
/*
 * The following four functions simply calculate the result while operator functions will
 * consider the accuracy.
 */
Numerical add(const Numerical& n1, const Numerical& n2) {
    if(n1.isZero())
        return Numerical(n2);
    else if(n2.isZero())
        return Numerical(n1);
    else if ((n1.length ^ n2.length) < 0) { // NOLINT(hicpp-signed-bitwise)
        if (n1.length > 0) {
            Numerical shallow_copy(n2.byte, -n2.length, n2.power);
            auto result = sub(n1, shallow_copy);
            shallow_copy.byte = nullptr;
            Q_UNUSED(shallow_copy)
            return result;
        }
        else {
            Numerical shallow_copy(n1.byte, -n1.length, n1.power);
            auto result = sub(n2, shallow_copy);
            shallow_copy.byte = nullptr;
            Q_UNUSED(shallow_copy)
            return result;
        }
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
        int length = (big->power + std::max(bigSize - big->power, smallSize - small->power));
        auto byte = (unsigned long*)malloc(length * sizeof(long));
        memcpy(byte + length - bigSize, big->byte, bigSize * sizeof(long));
        memset(byte, 0, (length - bigSize) * sizeof(long));

        int index = length - big->power + small->power - smallSize;
        unsigned long aByte;
        unsigned long carry = 0;
        unsigned long carry_temp;
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
            byte = (unsigned long*)realloc(byte, length * sizeof(long));
            byte[length - 1] = 1;
        }
        ////////////////////////////////////Out put////////////////////////////////////////
        if(big->length < 0)
            length = -length;
        return Numerical(byte, length, power);
    }
}

Numerical sub(const Numerical& n1, const Numerical& n2) {
    if(n1.isZero())
        return -n2;
    else if(n2.isZero())
        return Numerical(n1);
    else if (n1.length > 0) {
        if (n2.length < 0) {
            Numerical shallow_copy(n2.byte, -n2.length, n2.power);
            Numerical result = add(n1, shallow_copy);
            shallow_copy.byte = nullptr;
            Q_UNUSED(shallow_copy)
            return result;
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
            int length = (big->power + std::max(bigSize - big->power, smallSize - small->power));
            auto byte = (unsigned long*)malloc(length * sizeof(long));
            memcpy(byte + length - bigSize, big->byte, bigSize * sizeof(long));
            memset(byte, 0, (length - bigSize) * sizeof(long));

            int index = length - big->power + small->power - smallSize;
            unsigned long aByte;
            unsigned long carry = 0;
            unsigned long carry_temp;
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
                length = -length;
            Numerical result(byte, length, big->power);
            cutZero(result);
            return result;
        }
    }
    else {
        Numerical shallow_copy(n1.byte, -n1.length, n1.power);
        if (n2.length > 0) {
            Numerical result = add(shallow_copy, n2);
            result.toOpposite();
            shallow_copy.byte = nullptr;
            Q_UNUSED(shallow_copy)
            return result;
        }
        else {
            Numerical shallow_copy_1(n2.byte, n2.length, n2.power);
            Numerical result = sub(shallow_copy_1, shallow_copy);
            shallow_copy.byte = shallow_copy_1.byte = nullptr;
            Q_UNUSED(shallow_copy)
            Q_UNUSED(shallow_copy_1)
            return result;
        }
    }
}

Numerical mul(const Numerical& n1, const Numerical& n2) {
    if(n1 == basicConst->get_1())
        return Numerical(n2);
    else if(n2 == basicConst->get_1())
        return Numerical(n1);
    else {
        const int size1 = n1.getSize();
        const int size2 = n2.getSize();
        const int last2 = size2 - 1;
        //Estimate the ed of result first. we will calculate it accurately later.
        auto length = size1 + size2 - 1;
        auto byte = (unsigned long*)calloc(length, sizeof(long));
        unsigned long aByte;
        unsigned long carry = 0;
        //Every time the outer loop finished once, the position of carry will be reset. So we have to save the data.
        unsigned long carry_last = 0;
        unsigned long carry_temp;
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
            byte = (unsigned long*)realloc(byte, length * sizeof(long));
            byte[length - 1] = carry_last;
        }
        ////////////////////////////////////Out put////////////////////////////////////////
        if((n1.length ^ n2.length) < 0) // NOLINT(hicpp-signed-bitwise)
            length = -length;
        return Numerical(byte, length, power);
    }
}

Numerical div(const Numerical& n1, const Numerical& n2) {
    if(Q_UNLIKELY(n2.isZero()))
        qFatal("Divide by zero!");

    if(!n1.isZero()) {
        if(n2 != basicConst->get_1()) {
            Numerical n1_copy(n1);
            Numerical n2_copy(n2);
            n1_copy.toAbs();
            n2_copy.toAbs();
            n1_copy.a = n2_copy.a = 0;
            ////////////////////////////////Calculate cursory first//////////////////////////////////////
            //Estimate the length of result.
            int length = basicConst->getGlobalPrecision();
            //let n1_copy's power equal to n2_copy, power of the result will change correspondingly.
            int power = n1.getPower() - n2.getPower();
            n1_copy.power = n2.power;
            auto byte = (unsigned long*)calloc(length, sizeof(long));

            auto temp_arr = (unsigned long*)malloc(sizeof(long));
            Numerical temp(temp_arr, 1, 0);
            for (int i = length - 1; i >= 0; --i) {
                unsigned long large = ULONG_MAX;
                temp_arr[0] = ULONG_MAX / 2;
                unsigned long small = 0;
                while(true) {
                    Numerical mul = temp * n2_copy;
                    if(mul < n1_copy)
                        small = temp_arr[0];
                    else if (mul > n1_copy)
                        large = temp_arr[0];
                    else {
                        byte[i] = temp_arr[0];
                        goto stop;
                    }

                    temp_arr[0] = small + large;
                    if(temp_arr[0] < small) {
                        temp_arr[0] /= 2;
                        temp_arr[0] += (ULONG_MAX / 2 + 1);
                    }
                    else
                        temp_arr[0] /= 2;

                    if(small + 1 == large) {
                        n1_copy -= temp * n2_copy;
                        break;
                    }
                }
                ++n1_copy.power;
                byte[i] = temp_arr[0];
            }
            ////////////////////////////////////Out put////////////////////////////////////////
            stop:
            if((n1.length ^ n2.length) < 0) // NOLINT(hicpp-signed-bitwise)
                length = -length;
            //1 comes from the algorithm
            Numerical result(byte, length, power, 1);
            cutZero(result);
            return result;
        }
        else
            return Numerical(n1);
    }
    else
        return getZero();
}
/*
 * If the length of new array is larger than GlobalPrecision, it will be set to GlobalPrecision.
 * Return true if array is cut.
 */
bool cutLength(Numerical& n) {
    bool result = false;
    int size = n.getSize();

    if(size > basicConst->getGlobalPrecision()) {
        result = true;
        int cutFrom = size - basicConst->getGlobalPrecision();
        auto new_byte = (unsigned long*)malloc(basicConst->getGlobalPrecision() * sizeof(long));
        memcpy(new_byte, n.byte + cutFrom, basicConst->getGlobalPrecision() * sizeof(long));
        free(n.byte);
        n.byte = new_byte;
        auto length = basicConst->getGlobalPrecision();
        if(n.length < 0)
            length = -length;
        n.length = length;
        n.a = n.a != 0;
    }
    return result;
}
/*
 * Cut zeros from the beginning.
 */
void cutZero(Numerical& n) {
    int size = n.getSize();
    int id = size - 1;
    while(n.byte[id] == 0 && id > 0)
        --id;
    ++id;

    if(id != size) {
        int shorten = size - id;
        n.byte = (unsigned long*)realloc(n.byte, id * sizeof(long));
        size = id;
        if(n.length < 0)
            size = -size;
        n.length = size;

        if(n.byte[id - 1] != 0)
            n.power -= shorten;
        else
            n.power = 0;
    }
}