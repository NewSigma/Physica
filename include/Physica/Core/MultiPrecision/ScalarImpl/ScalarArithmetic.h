/*
 * Copyright 2020-2023 Weibo He.
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
#pragma once

#include "Physica/Core/MultiPrecision/BasicImpl/AddBasic.h"
#include "Physica/Core/MultiPrecision/BasicImpl/DivBasic.h"
#include "Physica/Core/MultiPrecision/BasicImpl/Util/ArraySupport.h"
#include "Physica/Core/MultiPrecision/BasicImpl/Util/Bitwise.h"
#include "Physica/Core/Exception/DivideByZeroException.h"

namespace Physica::Core {
    ///////////////////////////////////////BasicCalculates////////////////////////////////////////////
    inline Scalar<MultiPrecision> Scalar<MultiPrecision>::add(const Scalar<MultiPrecision>& s1, const Scalar<MultiPrecision>& s2) {
        if(s1.isZero())
            return Scalar<MultiPrecision>(static_cast<const Scalar<MultiPrecision>&>(s2));
        else if(s2.isZero())
            return Scalar<MultiPrecision>(static_cast<const Scalar<MultiPrecision>&>(s1));
        else if (!matchSign(s1, s2)) {
            if (s1.length > 0) {
                Scalar shallow_copy(const_cast<MPUnit*>(s2.byte), -s2.length, s2.power);
                auto result = sub(s1, shallow_copy);
                shallow_copy.byte = nullptr;
                return result;
            }
            else {
                Scalar shallow_copy(const_cast<MPUnit*>(s1.byte), -s1.length, s1.power);
                auto result = sub(s2, shallow_copy);
                shallow_copy.byte = nullptr;
                return result;
            }
        }
        else {
            const Scalar* big;
            const Scalar* small;
            if (s1.power > s2.power) {
                big = &s1;
                small = &s2;
            }
            else {
                big = &s2;
                small = &s1;
            }
            const int bigSize = big->getSize();
            const int smallSize = small->getSize();
            //Estimate the ed of result first, will calculate it accurately later.
            int length = (big->power + std::max(bigSize - big->power, smallSize - small->power));
            length = length > GlobalPrecision
                     ? GlobalPrecision : length;
            auto byte = reinterpret_cast<MPUnit*>(malloc(length * sizeof(MPUnit)));
            /* Init byte */ {
                const auto copySize = bigSize > length ? length : bigSize;
                const auto clearSize = length - copySize;
                memset(byte, 0, clearSize * sizeof(MPUnit));
                memcpy(byte + clearSize, big->byte + bigSize - copySize, copySize * sizeof(MPUnit));
            }
            bool carry;
            /* Add and carry */ {
                //usableSmall is the part whose add result will fall in GlobalPrecision.
                int usableSmall = small->power - (big->power - GlobalPrecision);
                MPUnit a = usableSmall < 0;
                usableSmall = usableSmall > smallSize
                               ? smallSize : (a ? 0 : usableSmall);
                carry = addArrWithArrEq(small->byte + smallSize - usableSmall
                        , byte + length + small->power - big->power - usableSmall, usableSmall);
                //usableSmall is also the index which we should carry to.
                while(carry != 0 && usableSmall < length) {
                    MPUnit temp = byte[usableSmall] + 1;
                    byte[usableSmall] = temp;
                    carry = temp < carry;
                    ++usableSmall;
                }
            }
            ///////////////////////////////////////Get byte, length and power//////////////////////////
            int power = big->power;
            if(carry) {
                ++length;
                ++power;
                byte = reinterpret_cast<MPUnit*>(realloc(byte, length * sizeof(MPUnit)));
                byte[length - 1] = 1;
            }
            ////////////////////////////////////Out put////////////////////////////////////////
            return Scalar<MultiPrecision>(byte, big->length < 0 ? -length : length, power);
        }
    }

    inline Scalar<MultiPrecision> Scalar<MultiPrecision>::sub(const Scalar<MultiPrecision>& s1, const Scalar<MultiPrecision>& s2) {
        if(s1.isZero())
            return Scalar<MultiPrecision>(static_cast<Scalar<MultiPrecision>&&>(-s2));
        else if(s2.isZero())
            return Scalar<MultiPrecision>(static_cast<const Scalar<MultiPrecision>&>(s1));
        else if (s1.length > 0) {
            if (s2.length < 0) {
                Scalar shallow_copy(const_cast<MPUnit*>(s2.byte), -s2.length, s2.power);
                Scalar result = add(s1, shallow_copy);
                shallow_copy.byte = nullptr;
                return result;
            }
            else {
                const Scalar* big;
                const Scalar* small;
                bool changeSign = false;
                if (s1.power > s2.power) {
                    big = &s1;
                    small = &s2;
                }
                else {
                    changeSign = true;
                    big = &s2;
                    small = &s1;
                }
                redo:
                const int bigSize = big->getSize();
                const int smallSize = small->getSize();
                //Estimate the ed of result first, will calculate it accurately later.
                int length = (big->power + std::max(bigSize - big->power, smallSize - small->power));
                length = length > GlobalPrecision
                         ? GlobalPrecision : length;
                auto byte = reinterpret_cast<MPUnit*>(malloc(length * sizeof(MPUnit)));
                /* Init byte */ {
                    const auto copySize = bigSize > length ? length : bigSize;
                    const auto clearSize = length - copySize;
                    memset(byte, 0, clearSize * sizeof(MPUnit));
                    memcpy(byte + clearSize, big->byte + bigSize - copySize, copySize * sizeof(MPUnit));
                }
                bool carry;
                MPUnit a;
                /* Sub and carry */ {
                    //usableSmall is the part whose sub result will fall in GlobalPrecision.
                    int usableSmall = small->power - (big->power - GlobalPrecision);
                    a = usableSmall < 0;
                    usableSmall = usableSmall > smallSize
                                   ? smallSize : (a ? 0 : usableSmall);
                    int smallEnd = length + small->power - big->power;
                    carry = subArrByArrEq(byte + smallEnd - usableSmall
                            , small->byte + smallSize - usableSmall, usableSmall);
                    //usableSmall is also the index which we should carry to.
                    MPUnit temp1, temp2;
                    while(carry != 0 && smallEnd < length) {
                        temp1 = byte[smallEnd];
                        temp2 = temp1 - 1;
                        byte[smallEnd] = temp2;
                        carry = temp1 < temp2;
                        ++smallEnd;
                    }
                }

                if(carry) {
                    auto temp = big;
                    big = small;
                    small = temp;
                    changeSign = !changeSign;
                    free(byte);
                    goto redo;
                }
                Scalar<MultiPrecision> result(byte, changeSign ? -length : length, big->power);
                result.cutZero();
                return result;
            }
        }
        else {
            Scalar shallow_copy(const_cast<MPUnit*>(s1.byte), -s1.length, s1.power);
            if (s2.length > 0) {
                Scalar result = add(shallow_copy, s2);
                result.toOpposite();
                shallow_copy.byte = nullptr;
                return result;
            }
            else {
                Scalar shallow_copy_1(const_cast<MPUnit*>(s2.byte), -s2.length, s2.power);
                Scalar result = sub(shallow_copy_1, shallow_copy);
                shallow_copy.byte = shallow_copy_1.byte = nullptr;
                return result;
            }
        }
    }
    //Optimize: length may be too long and it is unnecessary, cut it and consider the accuracy.
    inline Scalar<MultiPrecision> Scalar<MultiPrecision>::mul(const Scalar<MultiPrecision>& s1, const Scalar<MultiPrecision>& s2) {
        if (s1.isZero() || s2.isZero())
            return Scalar<MultiPrecision>(0);
        if (s1 == BasicConst::getInstance()._1)
            return Scalar<MultiPrecision>(static_cast<const Scalar<MultiPrecision>&>(s2));
        if (s2 == BasicConst::getInstance()._1)
            return Scalar<MultiPrecision>(static_cast<const Scalar<MultiPrecision>&>(s1));
        const int size1 = s1.getSize();
        const int size2 = s2.getSize();
        //Estimate the ed of result first. we will calculate it accurately later.
        auto length = size1 + size2;
        auto byte = reinterpret_cast<MPUnit*>(calloc(length, sizeof(MPUnit)));
        for (int i = 0; i < size1; ++i)
            byte[i + size2] = mulAddArrByWord(byte + i, s2.byte, size2, s1.byte[i]);
        ///////////////////////////////////////Get byte, length and power//////////////////////////;
        int power = s1.power + s2.power + 1;
        if (byte[length - 1] == 0) {
            --length;
            --power;
            byte = reinterpret_cast<MPUnit*>(realloc(byte, length * sizeof(MPUnit)));
        }
        ////////////////////////////////////Out put////////////////////////////////////////
        return Scalar<MultiPrecision>(byte, matchSign(s1, s2) ? length : -length, power);
    }

    inline Scalar<MultiPrecision> Scalar<MultiPrecision>::div(const Scalar<MultiPrecision>& s1, const Scalar<MultiPrecision>& s2) {
        if(s2.isZero())
            throw DivideByZeroException();

        if(!s1.isZero()) {
            if(s2 != BasicConst::getInstance()._1) {
                const auto s1_size = s1.getSize(), s2_size = s2.getSize();
                //Add one to arr1_length to avoid precision loss during right shift.
                auto arr1_len = std::max(s1_size, s2_size) + 1;
                auto s1_blank = arr1_len - s1_size;
                auto arr1 = new MPUnit[arr1_len];
                memcpy(arr1 + s1_blank, s1.byte, s1_size * sizeof(MPUnit));
                memset(arr1, 0, s1_blank * sizeof(MPUnit));
                //Size of arr2 is arranged 1 less than arr1.
                auto arr2_len = arr1_len - 1;
                auto s2_blank = arr2_len - s2_size;
                auto arr2 = new MPUnit[arr2_len];
                memcpy(arr2 + s2_blank, s2.byte, s2_size * sizeof(MPUnit));
                memset(arr2, 0, s2_blank * sizeof(MPUnit));
                /*
                 * We shift s1 and s2, making the less highest bit of s1 is set and the highest bit of s2 is set
                 * to meet the acquirement of the function divArrByFullArrWith1Word().
                 */
                const int s1_shift = static_cast<int>(countLeadingZeros(s1.byte[s1_size - 1])) - 1;
                if(s1_shift > 0)
                    byteLeftShiftEq(arr1, arr1_len, s1_shift);
                else
                    byteRightShiftEq(arr1, arr1_len, -s1_shift);
                const int s2_shift = static_cast<int>(countLeadingZeros(s2.byte[s2_size - 1]));
                byteLeftShiftEq(arr2, arr2_len, s2_shift);
                ////////////////////////////////Calculate cursory first//////////////////////////////////////
                //Estimate the length of result.
                int length = GlobalPrecision;
                auto byte = reinterpret_cast<MPUnit*>(malloc(length * sizeof(MPUnit)));

                for(int i = length - 1; i >= 0; --i) {
                    byte[i] = divArrByFullArrWith1Word(arr1, arr2, arr2_len);
                    arr1[arr2_len] -= mulSubArrByWord(arr1, arr2, arr2_len, byte[i]);
                    byteLeftShiftEq(arr1, arr1_len, MPUnitWidth);
                }
                delete[] arr1;
                delete[] arr2;
                ////////////////////////////////////Out put////////////////////////////////////////
                return Scalar<MultiPrecision>(byte, matchSign(s1, s2) ? length : -length
                        , s1.getPower() - s2.getPower() - 1) >> (s1_shift - s2_shift);
            }
            else
                return static_cast<const Scalar<MultiPrecision>&>(s1);
        }
        else
            return Scalar<MultiPrecision>(static_cast<SignedMPUnit>(0));
    }
    /**
     * If the length of new array is larger than GlobalPrecision, it will be set to GlobalPrecision.
     * Return true if array is cut.
     */
    inline bool Scalar<MultiPrecision>::cutLength(Scalar<MultiPrecision>& s) {
        bool result = false;
        int size = s.getSize();

        if(size > GlobalPrecision) {
            result = true;
            int cutFrom = size - GlobalPrecision;
            auto new_byte = reinterpret_cast<MPUnit*>(malloc(GlobalPrecision * sizeof(MPUnit)));
            memcpy(new_byte, s.byte + cutFrom, GlobalPrecision * sizeof(MPUnit));
            free(s.byte);
            s.byte = new_byte;
            s.length = s.length > 0 ? GlobalPrecision : -GlobalPrecision;
        }
        return result;
    }
}
