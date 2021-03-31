/*
 * Copyright 2020-2021 WeiBo He.
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
#ifndef PHYSICA_SCALARARITHMETIC_H
#define PHYSICA_SCALARARITHMETIC_H

#include "Physica/Core/MultiPrecision/BasicImpl/AddBasic.h"
#include "Physica/Core/MultiPrecision/BasicImpl/DivBasic.h"
#include "Physica/Core/MultiPrecision/BasicImpl/Util/ArraySupport.h"
#include "Physica/Core/MultiPrecision/BasicImpl/Util/Bitwise.h"
#include "Physica/Core/Exception/DivideByZeroException.h"
/*!
 * This file is part of implementations of \Scalar.
 * Do not include this header file, include Scalar.h instead.
 */
//Forward declaration
namespace Physica::Core {
    inline bool operator!=(const Internal::AbstractScalar<MultiPrecision>& s1, const Internal::AbstractScalar<MultiPrecision>& s2);
}

namespace Physica::Core::Internal {
    ///////////////////////////////////////BasicCalculates////////////////////////////////////////////
    /*!
     * WithError functions will calculate the error while NoError functions will not.
     * They have slightly difference between each other.
     */
    inline Scalar<MultiPrecision, true> AbstractScalar<MultiPrecision>::addWithError (
            const AbstractScalar<MultiPrecision>& s1, const AbstractScalar<MultiPrecision>& s2) {
        if(s1.isZero())
            return Scalar<MultiPrecision, true>(static_cast<const Scalar<MultiPrecision, false>&>(s2));
        else if(s2.isZero())
            return Scalar<MultiPrecision, true>(static_cast<const Scalar<MultiPrecision, false>&>(s1));
        else if (!matchSign(s1, s2)) {
            if (s1.length > 0) {
                AbstractScalar shallow_copy(const_cast<MPUnit*>(s2.byte), -s2.length, s2.power);
                auto result = subWithError(s1, shallow_copy);
                shallow_copy.byte = nullptr;
                Q_UNUSED(shallow_copy)
                return result;
            }
            else {
                AbstractScalar shallow_copy(const_cast<MPUnit*>(s1.byte), -s1.length, s1.power);
                auto result = subWithError(s2, shallow_copy);
                shallow_copy.byte = nullptr;
                Q_UNUSED(shallow_copy)
                return result;
            }
        }
        else {
            const AbstractScalar* big;
            const AbstractScalar* small;
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
            MPUnit a;
            /* Add and carry */ {
                //usableSmall is the part whose add result will fall in GlobalPrecision.
                int usableSmall = small->power - (big->power - GlobalPrecision);
                a = usableSmall < 0;
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
            return Scalar<MultiPrecision, true>(byte, big->length < 0 ? -length : length, power, a);
        }
    }

    inline Scalar<MultiPrecision, false> AbstractScalar<MultiPrecision>::addNoError (
            const AbstractScalar<MultiPrecision>& s1, const AbstractScalar<MultiPrecision>& s2) {
        if(s1.isZero())
            return Scalar<MultiPrecision, false>(static_cast<const Scalar<MultiPrecision, false>&>(s2));
        else if(s2.isZero())
            return Scalar<MultiPrecision, false>(static_cast<const Scalar<MultiPrecision, false>&>(s1));
        else if (!matchSign(s1, s2)) {
            if (s1.length > 0) {
                AbstractScalar shallow_copy(const_cast<MPUnit*>(s2.byte), -s2.length, s2.power);
                auto result = subNoError(s1, shallow_copy);
                shallow_copy.byte = nullptr;
                Q_UNUSED(shallow_copy)
                return result;
            }
            else {
                AbstractScalar shallow_copy(const_cast<MPUnit*>(s1.byte), -s1.length, s1.power);
                auto result = subNoError(s2, shallow_copy);
                shallow_copy.byte = nullptr;
                Q_UNUSED(shallow_copy)
                return result;
            }
        }
        else {
            const AbstractScalar* big;
            const AbstractScalar* small;
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
            return Scalar<MultiPrecision, false>(byte, big->length < 0 ? -length : length, power);
        }
    }

    inline Scalar<MultiPrecision, true> AbstractScalar<MultiPrecision>::subWithError(
            const AbstractScalar<MultiPrecision>& s1, const AbstractScalar<MultiPrecision>& s2) {
        if(s1.isZero())
            return Scalar<MultiPrecision, true>(static_cast<Scalar<MultiPrecision, false>&&>(-s2));
        else if(s2.isZero())
            return Scalar<MultiPrecision, true>(static_cast<const Scalar<MultiPrecision, false>&>(s1));
        else if (s1.length > 0) {
            if (s2.length < 0) {
                AbstractScalar shallow_copy(const_cast<MPUnit*>(s2.byte), -s2.length, s2.power);
                auto result = addWithError(s1, shallow_copy);
                shallow_copy.byte = nullptr;
                Q_UNUSED(shallow_copy)
                return result;
            }
            else {
                const AbstractScalar* big;
                const AbstractScalar* small;
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
                    //byte[smallEnd] is the element of \byte, which corresponds the end of small->byte.
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
                Scalar<MultiPrecision, true> result(byte, changeSign ? -length : length, big->power, a);
                cutZero(result);
                return result;
            }
        }
        else {
            AbstractScalar shallow_copy(const_cast<MPUnit*>(s1.byte), -s1.length, s1.power);
            if (s2.length > 0) {
                auto result = addWithError(shallow_copy, s2);
                result.toOpposite();
                shallow_copy.byte = nullptr;
                Q_UNUSED(shallow_copy)
                return result;
            }
            else {
                AbstractScalar shallow_copy_1(const_cast<MPUnit*>(s2.byte), -s2.length, s2.power);
                auto result = subWithError(shallow_copy_1, shallow_copy);
                shallow_copy.byte = shallow_copy_1.byte = nullptr;
                Q_UNUSED(shallow_copy)
                Q_UNUSED(shallow_copy_1)
                return result;
            }
        }
    }

    inline Scalar<MultiPrecision, false> AbstractScalar<MultiPrecision>::subNoError(
            const AbstractScalar<MultiPrecision>& s1, const AbstractScalar<MultiPrecision>& s2) {
        if(s1.isZero())
            return Scalar<MultiPrecision, false>(static_cast<Scalar<MultiPrecision, false>&&>(-s2));
        else if(s2.isZero())
            return Scalar<MultiPrecision, false>(static_cast<const Scalar<MultiPrecision, false>&>(s1));
        else if (s1.length > 0) {
            if (s2.length < 0) {
                AbstractScalar shallow_copy(const_cast<MPUnit*>(s2.byte), -s2.length, s2.power);
                Scalar result = addNoError(s1, shallow_copy);
                shallow_copy.byte = nullptr;
                Q_UNUSED(shallow_copy)
                return result;
            }
            else {
                const AbstractScalar* big;
                const AbstractScalar* small;
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
                Scalar<MultiPrecision, false> result(byte, changeSign ? -length : length, big->power);
                cutZero(result);
                return result;
            }
        }
        else {
            AbstractScalar shallow_copy(const_cast<MPUnit*>(s1.byte), -s1.length, s1.power);
            if (s2.length > 0) {
                Scalar result = addNoError(shallow_copy, s2);
                result.toOpposite();
                shallow_copy.byte = nullptr;
                Q_UNUSED(shallow_copy)
                return result;
            }
            else {
                AbstractScalar shallow_copy_1(const_cast<MPUnit*>(s2.byte), -s2.length, s2.power);
                Scalar result = subNoError(shallow_copy_1, shallow_copy);
                shallow_copy.byte = shallow_copy_1.byte = nullptr;
                Q_UNUSED(shallow_copy)
                Q_UNUSED(shallow_copy_1)
                return result;
            }
        }
    }

    inline Scalar<MultiPrecision, true> AbstractScalar<MultiPrecision>::mulWithError(
            const AbstractScalar<MultiPrecision>& s1, const AbstractScalar<MultiPrecision>& s2) {
        if (s1.isZero() || s2.isZero())
            return Scalar<MultiPrecision, true>(0);
        if (s1 == BasicConst::getInstance()._1)
            return Scalar<MultiPrecision, true>(static_cast<const Scalar<MultiPrecision, false>&>(s2));
        if (s2 == BasicConst::getInstance()._1)
            return Scalar<MultiPrecision, true>(static_cast<const Scalar<MultiPrecision, false>&>(s1));
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
        return Scalar<MultiPrecision, true>(byte, matchSign(s1, s2) ? length : -length, power);
    }
    //Optimize: length may be too long and it is unnecessary, cut it and consider the accuracy.
    inline Scalar<MultiPrecision, false> AbstractScalar<MultiPrecision>::mulNoError(
            const AbstractScalar<MultiPrecision>& s1, const AbstractScalar<MultiPrecision>& s2) {
        if (s1.isZero() || s2.isZero())
            return Scalar<MultiPrecision, true>(0);
        if (s1 == BasicConst::getInstance()._1)
            return Scalar<MultiPrecision, false>(static_cast<const Scalar<MultiPrecision, false>&>(s2));
        if (s2 == BasicConst::getInstance()._1)
            return Scalar<MultiPrecision, false>(static_cast<const Scalar<MultiPrecision, false>&>(s1));
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
        return Scalar<MultiPrecision, false>(byte, matchSign(s1, s2) ? length : -length, power);
    }

    inline Scalar<MultiPrecision, true> AbstractScalar<MultiPrecision>::divWithError(
            const AbstractScalar<MultiPrecision>& s1, const AbstractScalar<MultiPrecision>& s2) {
        if(Q_UNLIKELY(s2.isZero()))
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
                //Accuracy 1 comes from the algorithm
                return Scalar<MultiPrecision, true>(byte, matchSign(s1, s2) ? length : -length
                        , s1.getPower() - s2.getPower() - 1, 1) >> (s1_shift - s2_shift);
            }
            else
                return Scalar<MultiPrecision, true>(static_cast<const Scalar<MultiPrecision, false>&>(s1));
        }
        else
            return Scalar<MultiPrecision, true>(static_cast<SignedMPUnit>(0));
    }

    inline Scalar<MultiPrecision, false> AbstractScalar<MultiPrecision>::divNoError(
            const AbstractScalar<MultiPrecision>& s1, const AbstractScalar<MultiPrecision>& s2) {
        if(Q_UNLIKELY(s2.isZero()))
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
                return Scalar<MultiPrecision, false>(byte, matchSign(s1, s2) ? length : -length
                        , s1.getPower() - s2.getPower() - 1) >> (s1_shift - s2_shift);
            }
            else
                return static_cast<const Scalar<MultiPrecision, false>&>(s1);
        }
        else
            return Scalar<MultiPrecision, false>(static_cast<SignedMPUnit>(0));
    }
    /*!
     * If the length of new array is larger than GlobalPrecision, it will be set to GlobalPrecision.
     * Return true if array is cut.
     */
    template<>
    inline bool AbstractScalar<MultiPrecision>::cutLength(Scalar<MultiPrecision, true>& s) {
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
            s.a = s.a != 0;
        }
        return result;
    }

    template<>
    inline bool AbstractScalar<MultiPrecision>::cutLength(Scalar<MultiPrecision, false>& s) {
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

#endif
