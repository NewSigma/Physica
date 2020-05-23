/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_ARRAYSUPPORT_H
#define PHYSICA_ARRAYSUPPORT_H

#include <cstring>
#include <memory>
#include "Core/Header/Numerical.h"
#include "Bitwise.h"

namespace Physica::Core {
    // operator<<
    [[nodiscard]] inline std::unique_ptr<NumericalUnit[]> byteLeftShift(const NumericalUnit* byte, unsigned int length, unsigned int shift) {
        const unsigned int quotient = shift / NumericalUnitWidth;
        auto result = std::unique_ptr<NumericalUnit[]>(new NumericalUnit[length]);
        if(quotient < length) {
            if(quotient != 0) {
                memcpy(result.get() + quotient, byte, (length - quotient) * sizeof(NumericalUnit));
                memset(result.get(), 0, quotient * sizeof(NumericalUnit));
                shift -= quotient * NumericalUnitWidth;
            }
            if(shift != 0) {
                NumericalUnit carry = 0, temp;
                for(unsigned int i = quotient; i < length - 1; ++i) {
                    temp = result[i] >> (NumericalUnitWidth - shift);
                    result[i] <<= shift;
                    result[i] |= carry;
                    carry = temp;
                }
                result[length - 1] <<= shift;
                result[length - 1] |= carry;
            }
        }
        return result;
    }
    // operator>>
    [[nodiscard]] inline std::unique_ptr<NumericalUnit[]> byteRightShift(const NumericalUnit* byte, size_t length, size_t shift) {
        const unsigned int quotient = shift / NumericalUnitWidth;
        auto result = std::unique_ptr<NumericalUnit[]>(new NumericalUnit[length]);
        if(quotient < length) {
            auto bufferSize = length - quotient;
            if(quotient != 0) {
                memcpy(result.get(), byte + quotient, (length - quotient) * sizeof(NumericalUnit));
                memset(result.get() + bufferSize, 0, quotient * sizeof(NumericalUnit));
                shift -= quotient * NumericalUnitWidth;
            }
            if(shift != 0) {
                NumericalUnit carry = 0, temp;
                for(unsigned int i = bufferSize - 1; i > 0; --i) {
                    temp = result[i] << (NumericalUnitWidth - shift);
                    result[i] >>= shift;
                    result[i] |= carry;
                    carry = temp;
                }
                result[0] >>= shift;
                result[0] |= carry;
            }
        }
        return result;
    }
    // operator<<=
    inline void byteLeftShiftEq(NumericalUnit* byte, unsigned int length, unsigned int shift) {
        const unsigned int quotient = shift / NumericalUnitWidth;
        if(quotient < length) {
            if(quotient != 0) {
                memmove(byte + quotient, byte, (length - quotient) * sizeof(NumericalUnit));
                memset(byte, 0, quotient * sizeof(NumericalUnit));
                shift -= quotient * NumericalUnitWidth;
            }
            if(shift != 0) {
                NumericalUnit carry = 0, temp;
                for(unsigned int i = quotient; i < length - 1; ++i) {
                    temp = byte[i] >> (NumericalUnitWidth - shift);
                    byte[i] <<= shift;
                    byte[i] |= carry;
                    carry = temp;
                }
                byte[length - 1] <<= shift;
                byte[length - 1] |= carry;
            }
        }
    }
    // operator>>=
    inline void byteRightShiftEq(NumericalUnit* byte, size_t length, size_t shift) {
        const unsigned int quotient = shift / NumericalUnitWidth;
        if(quotient < length) {
            auto bufferSize = length - quotient;
            if(quotient != 0) {
                memmove(byte, byte + quotient, bufferSize * sizeof(NumericalUnit));
                memset(byte + bufferSize, 0, quotient * sizeof(NumericalUnit));
                shift -= quotient * NumericalUnitWidth;
            }
            if(shift != 0) {
                NumericalUnit carry = 0, temp;
                for(unsigned int i = bufferSize - 1; i > 0; --i) {
                    temp = byte[i] << (NumericalUnitWidth - shift);
                    byte[i] >>= shift;
                    byte[i] |= carry;
                    carry = temp;
                }
                byte[0] >>= shift;
                byte[0] |= carry;
            }
        }
    }
}

#endif