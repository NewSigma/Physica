/*
 * Copyright 2020 WeiBo He.
 *
 * This file is part of Physica.
 *
 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
 */
#ifndef PHYSICA_ARRAYSUPPORT_H
#define PHYSICA_ARRAYSUPPORT_H

#include <cstring>
#include <memory>

namespace Physica::Core {
    // operator<<
    [[nodiscard]] inline std::unique_ptr<MPUnit[]> byteLeftShift(
                        const MPUnit* __restrict byte, unsigned int length, unsigned int shift) {
        const unsigned int quotient = shift / MPUnitWidth;
        auto result = std::unique_ptr<MPUnit[]>(new MPUnit[length]);
        if(quotient < length) {
            if(quotient != 0) {
                memcpy(result.get() + quotient, byte, (length - quotient) * sizeof(MPUnit));
                memset(result.get(), 0, quotient * sizeof(MPUnit));
                shift -= quotient * MPUnitWidth;
            }
            if(shift != 0) {
                MPUnit carry = 0, temp;
                for(unsigned int i = quotient; i < length - 1; ++i) {
                    temp = result[i] >> (MPUnitWidth - shift);
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
    [[nodiscard]] inline std::unique_ptr<MPUnit[]> byteRightShift(
                            const MPUnit* __restrict byte, size_t length, size_t shift) {
        const unsigned int quotient = shift / MPUnitWidth;
        auto result = std::unique_ptr<MPUnit[]>(new MPUnit[length]);
        if(quotient < length) {
            auto bufferSize = length - quotient;
            if(quotient != 0) {
                memcpy(result.get(), byte + quotient, (length - quotient) * sizeof(MPUnit));
                memset(result.get() + bufferSize, 0, quotient * sizeof(MPUnit));
                shift -= quotient * MPUnitWidth;
            }
            if(shift != 0) {
                MPUnit carry = 0, temp;
                for(unsigned int i = bufferSize - 1; i > 0; --i) {
                    temp = result[i] << (MPUnitWidth - shift);
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
    inline void byteLeftShiftEq(MPUnit* __restrict byte, unsigned int length, unsigned int shift) {
        const unsigned int quotient = shift / MPUnitWidth;
        if(quotient < length) {
            if(quotient != 0) {
                memmove(byte + quotient, byte, (length - quotient) * sizeof(MPUnit));
                memset(byte, 0, quotient * sizeof(MPUnit));
                shift -= quotient * MPUnitWidth;
            }
            if(shift != 0) {
                MPUnit carry = 0, temp;
                for(unsigned int i = quotient; i < length - 1; ++i) {
                    temp = byte[i] >> (MPUnitWidth - shift);
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
    inline void byteRightShiftEq(MPUnit* __restrict byte, size_t length, size_t shift) {
        const unsigned int quotient = shift / MPUnitWidth;
        if(quotient < length) {
            auto bufferSize = length - quotient;
            if(quotient != 0) {
                memmove(byte, byte + quotient, bufferSize * sizeof(MPUnit));
                memset(byte + bufferSize, 0, quotient * sizeof(MPUnit));
                shift -= quotient * MPUnitWidth;
            }
            if(shift != 0) {
                MPUnit carry = 0, temp;
                for(unsigned int i = bufferSize - 1; i > 0; --i) {
                    temp = byte[i] << (MPUnitWidth - shift);
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