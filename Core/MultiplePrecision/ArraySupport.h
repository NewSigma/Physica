/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_ARRAYSUPPORT_H
#define PHYSICA_ARRAYSUPPORT_H

#include <cstring>
#include "Core/Header/Numerical.h"
#include "Bitwise.h"

namespace Physica::Core {
    inline void byteLeftShift(NumericalUnit* byte, unsigned int length, unsigned int shift) {
        const unsigned int quotient = shift / NumericalUnitWidth;
        bool overflow = quotient >= length;

        auto buffer = new NumericalUnit[length]{0};
        memcpy(buffer + quotient, byte, (overflow ? 0U : (length - quotient)) * sizeof(NumericalUnit));
        memcpy(byte, buffer, length * sizeof(NumericalUnit));
        delete[] buffer;

        if(!overflow) {
            const unsigned int remainder = shift - quotient * NumericalUnitWidth;
            unsigned int carry = 0, temp;
            for(unsigned int i = quotient; i < length - 1; ++i) {
                temp = byte[i] >> (NumericalUnitWidth - remainder);
                byte[i] <<= remainder;
                byte[i] |= carry;
                carry = temp;
            }
            byte[length - 1] <<= remainder;
            byte[length - 1] |= carry;
        }
    }

    inline void byteRightShift(NumericalUnit* byte, size_t length, size_t shift) {
        const unsigned int quotient = shift / NumericalUnitWidth;
        bool overflow = quotient >= length;

        auto buffer = new NumericalUnit[length]{0};
        memcpy(buffer, byte + quotient, (overflow ? 0U : (length - quotient)) * sizeof(NumericalUnit));
        memcpy(byte, buffer, length * sizeof(NumericalUnit));
        delete[] buffer;

        if(!overflow) {
            const unsigned int remainder = shift - quotient * NumericalUnitWidth;
            unsigned int carry = 0, temp;
            for(unsigned int i = length - quotient; i > 0; --i) {
                temp = byte[i] << (NumericalUnitWidth - remainder);
                byte[i] >>= remainder;
                byte[i] |= carry;
                carry = temp;
            }
            byte[length - 1] >>= remainder;
            byte[length - 1] |= carry;
        }
    }
}

#endif