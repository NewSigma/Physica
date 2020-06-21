/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_BASICCALC_H
#define PHYSICA_BASICCALC_H

#ifndef PHYSICA_SCALAR_H
    #include "Core/Header/MultiPrecition/Scalar.h"
#endif

namespace Physica::Core {
    ///////////////////////////////////////BasicCalculates////////////////////////////////////////////
    //Refactor: move the key algorithm into a separated method. We may not use every words.
    template<size_t maxPrecision>
    template<size_t maxPrecision2>
    inline Scalar<META_MAX<maxPrecision, maxPrecision2>::value, false> Scalar<maxPrecision, false>::add (
            const Scalar<maxPrecision, false>& s1, const Scalar<maxPrecision2, false>& s2) {
        static_assert(maxPrecision > 1 && maxPrecision2 > 1, "This function can only apply to multi-precision Scalar.");
        if(s1.isZero())
            return Scalar(s2);
        else if(s2.isZero())
            return Scalar(s1);
        else if ((s1.length ^ s2.length) < 0) { // NOLINT(hicpp-signed-bitwise)
            if (s1.length > 0) {
                Scalar shallow_copy(const_cast<ScalarUnit*>(s2.byte), -s2.length, s2.power);
                auto result = sub(s1, shallow_copy);
                shallow_copy.byte = nullptr;
                Q_UNUSED(shallow_copy)
                return result;
            }
            else {
                Scalar shallow_copy(const_cast<ScalarUnit*>(s1.byte), -s1.length, s1.power);
                auto result = sub(s2, shallow_copy);
                shallow_copy.byte = nullptr;
                Q_UNUSED(shallow_copy)
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
            length = length > maxPrecision
                     ? maxPrecision : length;
            auto byte = reinterpret_cast<ScalarUnit*>(malloc(length * sizeof(ScalarUnit)));
            /* Init byte */ {
                const auto copySize = bigSize > length ? length : bigSize;
                const auto clearSize = length - copySize;
                memset(byte, 0, clearSize * sizeof(ScalarUnit));
                memcpy(byte + clearSize, big->byte + bigSize - copySize, copySize * sizeof(ScalarUnit));
            }
            bool carry;
            ScalarUnit a;
            /* Add and carry */ {
                //useableSmall is the part whose add result will fall in maxPrecision.
                int useableSmall = small->power - (big->power - maxPrecision);
                a = useableSmall < 0;
                useableSmall = useableSmall > smallSize
                               ? smallSize : (a ? 0 : useableSmall);
                carry = addArrWithArrEq(small->byte + smallSize - useableSmall
                        , byte + length + small->power - big->power - useableSmall, useableSmall);
                //useableSmall is also the index which we should carry to.
                while(carry != 0 && useableSmall < length) {
                    ScalarUnit temp = byte[useableSmall] + 1;
                    byte[useableSmall] = temp;
                    carry = temp < carry;
                    ++useableSmall;
                }
            }
            ///////////////////////////////////////Get byte, length and power//////////////////////////
            int power = big->power;
            if(carry) {
                ++length;
                ++power;
                byte = reinterpret_cast<ScalarUnit*>(realloc(byte, length * sizeof(ScalarUnit)));
                byte[length - 1] = 1;
            }
            ////////////////////////////////////Out put////////////////////////////////////////
            return Scalar(byte, big->length < 0 ? -length : length, power, a);
        }
    }
    //Refactor: move the key algorithm into a separated method. We may not use every words.
    template<size_t maxPrecision>
    template<size_t maxPrecision2>
    inline Scalar<META_MAX<maxPrecision, maxPrecision2>::value, false> Scalar<maxPrecision, false>::sub(
            const Scalar<maxPrecision, false>& s1, const Scalar<maxPrecision2, false>& s2) {
        static_assert(maxPrecision > 1 && maxPrecision2 > 1, "This function can only apply to multi-precision Scalar.");
        if(s1.isZero())
            return -s2;
        else if(s2.isZero())
            return Scalar(s1);
        else if (s1.length > 0) {
            if (s2.length < 0) {
                Scalar shallow_copy(const_cast<ScalarUnit*>(s2.byte), -s2.length, s2.power);
                Scalar result = add(s1, shallow_copy);
                shallow_copy.byte = nullptr;
                Q_UNUSED(shallow_copy)
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
                length = length > maxPrecision
                         ? maxPrecision : length;
                auto byte = reinterpret_cast<ScalarUnit*>(malloc(length * sizeof(ScalarUnit)));
                /* Init byte */ {
                    const auto copySize = bigSize > length ? length : bigSize;
                    const auto clearSize = length - copySize;
                    memset(byte, 0, clearSize * sizeof(ScalarUnit));
                    memcpy(byte + clearSize, big->byte + bigSize - copySize, copySize * sizeof(ScalarUnit));
                }
                bool carry;
                ScalarUnit a;
                /* Add and carry */ {
                    //useableSmall is the part whose sub result will fall in maxPrecision.
                    int useableSmall = small->power - (big->power - maxPrecision);
                    a = useableSmall < 0;
                    useableSmall = useableSmall > smallSize
                                   ? smallSize : (a ? 0 : useableSmall);
                    carry = subArrByArrEq(byte + length + small->power - big->power - useableSmall
                            , small->byte + smallSize - useableSmall, useableSmall);
                    //useableSmall is also the index which we should carry to.
                    ScalarUnit temp1, temp2;
                    while(carry != 0 && useableSmall < length) {
                        temp1 = byte[useableSmall];
                        temp2 = temp1 - 1;
                        byte[useableSmall] = temp2;
                        carry = temp1 < temp2;
                        ++useableSmall;
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
                Scalar result(byte, changeSign ? -length : length, big->power, a);
                cutZero(result);
                return result;
            }
        }
        else {
            Scalar shallow_copy(const_cast<ScalarUnit*>(s1.byte), -s1.length, s1.power);
            if (s2.length > 0) {
                Scalar result = add(shallow_copy, s2);
                result.toOpposite();
                shallow_copy.byte = nullptr;
                Q_UNUSED(shallow_copy)
                return result;
            }
            else {
                Scalar shallow_copy_1(const_cast<ScalarUnit*>(s2.byte), s2.length, s2.power);
                Scalar result = sub(shallow_copy_1, shallow_copy);
                shallow_copy.byte = shallow_copy_1.byte = nullptr;
                Q_UNUSED(shallow_copy)
                Q_UNUSED(shallow_copy_1)
                return result;
            }
        }
    }

    template<size_t maxPrecision>
    template<size_t maxPrecision2>
    inline Scalar<META_MAX<maxPrecision, maxPrecision2>::value, false> Scalar<maxPrecision, false>::mul(
            const Scalar<maxPrecision, false>& s1, const Scalar<maxPrecision2, false>& s2) {
        static_assert(maxPrecision > 1 && maxPrecision2 > 1, "This function can only apply to multi-precision Scalar.");
        if(s1 == BasicConst::getInstance().get_1())
            return Scalar(s2);
        else if(s2 == BasicConst::getInstance().get_1())
            return Scalar(s1);
        else {
            const int size1 = s1.getSize();
            const int size2 = s2.getSize();
            //Estimate the ed of result first. we will calculate it accurately later.
            auto length = size1 + size2;
            auto byte = reinterpret_cast<ScalarUnit*>(calloc(length, sizeof(ScalarUnit)));
            for (int i = 0; i < size1; ++i)
                byte[i + size2] = mulAddArrByWord(byte + i, s2.byte, size2, s1.byte[i]);
            ///////////////////////////////////////Get byte, length and power//////////////////////////;
            int power = s1.power + s2.power + 1;
            if (byte[length - 1] == 0) {
                --length;
                --power;
                byte = reinterpret_cast<ScalarUnit*>(realloc(byte, length * sizeof(ScalarUnit)));
            }
            ////////////////////////////////////Out put////////////////////////////////////////
            return Scalar(byte, (s1.length ^ s2.length) < 0 ? -length : length, power); //NOLINT
        }
    }

    template<size_t maxPrecision>
    template<size_t maxPrecision2>
    inline Scalar<META_MAX<maxPrecision, maxPrecision2>::value, false> Scalar<maxPrecision, false>::div(
            const Scalar<maxPrecision, false>& s1, const Scalar<maxPrecision2, false>& s2) {
        static_assert(maxPrecision > 1 && maxPrecision2 > 1, "This function can only apply to multi-precision Scalar.");
        if(Q_UNLIKELY(s2.isZero()))
            qFatal("Encountered dividing by zero exception.");

        if(!s1.isZero()) {
            if(s2 != BasicConst::getInstance().get_1()) {
                const auto s1_size = s1.getSize(), s2_size = s2.getSize();
                //Add one to avoid precision loss during right shift.
                auto arr1_len = std::max(s1_size, s2_size) + 1;
                auto s1_blank = arr1_len - s1_size;
                auto arr1 = new ScalarUnit[arr1_len];
                memcpy(arr1 + s1_blank, s1.byte, s1_size * sizeof(ScalarUnit));
                memset(arr1, 0, s1_blank * sizeof(ScalarUnit));
                //Size of arr2 is arranged 1 less than arr1.
                auto arr2_len = arr1_len - 1;
                auto s2_blank = arr2_len - s2_size;
                auto arr2 = new ScalarUnit[arr2_len];
                memcpy(arr2 + s2_blank, s2.byte, s2_size * sizeof(ScalarUnit));
                memset(arr2, 0, s2_blank * sizeof(ScalarUnit));

                const int s1_shift = static_cast<int>(countLeadingZeros(s1.byte[s1_size - 1])) - 1;
                if(s1_shift > 0)
                    byteLeftShiftEq(arr1, arr1_len, s1_shift);
                else
                    byteRightShiftEq(arr1, arr1_len, -s1_shift);
                const int s2_shift = static_cast<int>(countLeadingZeros(s2.byte[s2_size - 1]));
                byteLeftShiftEq(arr2, arr2_len, s2_shift);
                ////////////////////////////////Calculate cursory first//////////////////////////////////////
                //Estimate the length of result.
                int length = maxPrecision;
                //let s1_copy's power equal to s2_copy, power of the result will change correspondingly.
                auto byte = reinterpret_cast<ScalarUnit*>(malloc(length * sizeof(ScalarUnit)));

                for(int i = length - 1; i >= 0; --i) {
                    byte[i] = divArrByFullArrWith1Word(arr1, arr2, arr2_len);
                    arr1[arr2_len] -= mulSubArrByWord(arr1, arr2, arr2_len, byte[i]);
                    byteLeftShiftEq(arr1, arr1_len, ScalarUnitWidth);
                }
                delete[] arr1;
                delete[] arr2;
                ////////////////////////////////////Out put////////////////////////////////////////
                //1 comes from the algorithm
                return Scalar(byte, (s1.length ^ s2.length) < 0 ? -length : length //NOLINT
                        , s1.getPower() - s2.getPower() - 1, 1) >> (s1_shift - s2_shift);
            }
            else
                return Scalar(s1);
        }
        else
            return Scalar(static_cast<SignedScalarUnit>(0));
    }
    /*!
     * If the length of new array is larger than maxPrecision, it will be set to maxPrecision.
     * Return true if array is cut.
     */
    template<size_t maxPrecision>
    inline bool Scalar<maxPrecision, false>::cutLength(Scalar<maxPrecision, false>& s) {
        static_assert(maxPrecision > 1, "This function can only apply to multi-precision Scalar.");
        bool result = false;
        int size = s.getSize();

        if(size > maxPrecision) {
            result = true;
            int cutFrom = size - maxPrecision;
            auto new_byte = reinterpret_cast<ScalarUnit*>(malloc(maxPrecision * sizeof(ScalarUnit)));
            memcpy(new_byte, s.byte + cutFrom, maxPrecision * sizeof(ScalarUnit));
            free(s.byte);
            s.byte = new_byte;
            s.length = s.length > 0 ? maxPrecision : -maxPrecision;
            s.a = s.a != 0;
        }
        return result;
    }
    /*!
     * Cut zeros from the beginning.
     */
    template<size_t maxPrecision>
    inline void Scalar<maxPrecision, false>::cutZero(Scalar<maxPrecision, false>& s) {
        static_assert(maxPrecision > 1, "This function can only apply to multi-precision Scalar.");
        const int size = s.getSize();
        int id = size - 1;
        while(s.byte[id] == 0 && id > 0)
            --id;
        ++id;

        if(id != size) {
            int shorten = size - id;
            s.byte = reinterpret_cast<ScalarUnit*>(realloc(s.byte, id * sizeof(ScalarUnit)));
            s.length = s.length > 0 ? id : -id;
            auto temp = s.power;
            s.power = s.byte[id - 1] != 0 ? (temp - shorten) : 0;
        }
    }
}

#endif