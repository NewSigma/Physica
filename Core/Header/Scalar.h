/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_SCALAR_H
#define PHYSICA_SCALAR_H

#include <iostream>
#include <cstring>
#include <QtCore/qlogging.h>
#include "Core/Header/SystemBits.h"
#include "Core/Header/Const.h"
#include "Core/MultiplePrecision/Basic/AddBasic.h"
#include "Core/MultiplePrecision/Basic/DivBasic.h"
#include "Core/MultiplePrecision/Util/ArraySupport.h"
#include "Core/Header/ElementaryFunction.h"

namespace Physica::Core {
    enum ScalarType {
        MultiPrecition,
        Float,
        Double
    };

    class Scalar {
        //Store effective digits.
        ScalarUnit* __restrict byte;
        /*
         * Length of byte = abs(length).
         * sign of length and sign of Scalar are same. (when Scalar != 0)
        */
        int length;
        /*
         * Number = (x0 +- a * (2 ^ __WORDSIZE) ^ (1-length)) * (2 ^ __WORDSIZE) ^power
         * We have not considered overflow of power in our codes elsewhere.
         */
        int power;
        //Accuracy
        ScalarUnit a;
    public:
        Scalar() noexcept;
        Scalar(int length, int power, ScalarUnit a = 0) noexcept;
        Scalar(const Scalar& n) noexcept;
        Scalar(Scalar&& n) noexcept;
        explicit Scalar(SignedScalarUnit unit, ScalarUnit a = 0) noexcept;
        explicit Scalar(double d, ScalarUnit a = 0);
        explicit Scalar(const char* s, ScalarUnit a = 0);
        explicit Scalar(const wchar_t* s, ScalarUnit a = 0);
        explicit Scalar(const std::string& s, ScalarUnit a = 0);
        explicit Scalar(const std::wstring& s, ScalarUnit a = 0);
        ~Scalar();
        /* Operators */
        explicit operator double() const;
        friend std::ostream& operator<<(std::ostream& os, const Scalar& n);
        Scalar operator+(const Scalar& n) const;
        Scalar operator-(const Scalar& n) const;
        Scalar operator*(const Scalar& n) const;
        Scalar operator/(const Scalar& n) const;
        Scalar operator<<(int bits) const;
        Scalar operator>>(int bits) const;
        void operator<<=(int bits) noexcept { *this = *this << bits; }
        void operator>>=(int bits) noexcept { *this = *this >> bits; }
        ScalarUnit& operator[](unsigned int index) { return byte[index]; }
        const ScalarUnit& operator[](unsigned int index) const { return byte[index]; }
        Scalar& operator=(const Scalar& n);
        Scalar& operator=(Scalar&& n) noexcept;
        Scalar operator^(const Scalar& n) const;
        Scalar operator-() const;
        /* Helpers */
        Scalar& applyError(const Scalar& error);
        void swap(Scalar& n) noexcept;
        Scalar& toAbs() noexcept { length = getSize(); return *this; }
        Scalar& toOpposite() noexcept { length = -length; return *this; }
        Scalar& toUnitA() noexcept { a = 1; return *this; }
        Scalar& clearA() noexcept { a = 0; return *this; }
        /* Getters */
        [[nodiscard]] int getLength() const noexcept { return length; }
        [[nodiscard]] int getPower() const noexcept { return power; }
        [[nodiscard]] ScalarUnit getA() const noexcept { return a; }
        [[nodiscard]] int getSize() const noexcept { return abs(length); }
        [[nodiscard]] bool isZero() const { return byte[getSize() - 1] == 0; }
        [[nodiscard]] bool isPositive() const { return !isZero() && length > 0; }
        [[nodiscard]] bool isNegative() const { return !isZero() && length < 0; }
        [[nodiscard]] bool isInteger() const { return getSize() == power + 1; }
        [[nodiscard]] Scalar getAccuracy() const;
        [[nodiscard]] Scalar getMaximum() const { return Scalar::add(*this, getAccuracy()).clearA(); }
        [[nodiscard]] Scalar getMinimum() const { return Scalar::sub(*this, getAccuracy()).clearA(); }
    protected:
        Scalar(ScalarUnit* byte, int length, int power, ScalarUnit a = 0);
        inline static Scalar add (const Scalar& n1, const Scalar& n2);
        inline static Scalar sub (const Scalar& n1, const Scalar& n2);
        inline static Scalar mul (const Scalar& n1, const Scalar& n2);
        inline static Scalar div (const Scalar& n1, const Scalar& n2);
        inline static bool cutLength(Scalar& n);
        inline static void cutZero(Scalar& n);
        /* Friends */
        friend class Solve;
        friend Scalar square(const Scalar& n);
        friend Scalar sqrt_light(const Scalar& n);
        friend Scalar ln_light(const Scalar& n);
    };
    bool absCompare(const Scalar& n1, const Scalar& n2);
    /* Operators */
    bool operator>(const Scalar& n1, const Scalar& n2);
    bool operator<(const Scalar& n1, const Scalar& n2);
    bool operator>=(const Scalar& n1, const Scalar& n2);
    bool operator<=(const Scalar& n1, const Scalar& n2);
    bool operator==(const Scalar& n1, const Scalar& n2);
    bool operator!=(const Scalar& n1, const Scalar& n2);
    /* Inline Implementations */
    inline Scalar operator+(const Scalar& n) { return Scalar(n); }
    inline void operator+=(Scalar& n1, const Scalar& n2) { n1 = n1 + n2; }
    inline void operator-=(Scalar& n1, const Scalar& n2) { n1 = n1 - n2; }
    inline void operator*=(Scalar& n1, const Scalar& n2) { n1 = n1 * n2; }
    inline void operator/=(Scalar& n1, const Scalar& n2) { n1 = n1 / n2; }
    inline void operator^=(Scalar& n1, const Scalar& n2) { n1 = n1 ^ n2; }
    inline Scalar& operator++(Scalar& n) { n += BasicConst::getInstance().get_1(); return n; }
    inline Scalar& operator--(Scalar& n) { n -= BasicConst::getInstance().get_1(); return n; }
    inline Scalar operator++(Scalar& n, int) { Scalar temp(n); n += BasicConst::getInstance().get_1(); return temp; } //NOLINT
    inline Scalar operator--(Scalar& n, int) { Scalar temp(n); n -= BasicConst::getInstance().get_1(); return temp; } //NOLINT
    inline Scalar getZero() { return Scalar(BasicConst::getInstance().get_0()); }
    inline Scalar getOne() { return Scalar(BasicConst::getInstance().get_1()); }
    inline Scalar getTwo() { return Scalar(BasicConst::getInstance().get_2()); }
    inline void swap(Scalar& n1, Scalar& n2) noexcept { n1.swap(n2); }
    /*
     * The following four functions simply calculate the result while operator functions will
     * consider the accuracy.
     */
    //Refactor: move the key algorithm into a separated method. We may not use every words.
    inline Scalar Scalar::add(const Scalar& n1, const Scalar& n2) {
        if(n1.isZero())
            return Scalar(n2);
        else if(n2.isZero())
            return Scalar(n1);
        else if ((n1.length ^ n2.length) < 0) { // NOLINT(hicpp-signed-bitwise)
            if (n1.length > 0) {
                Scalar shallow_copy(const_cast<ScalarUnit*>(n2.byte), -n2.length, n2.power);
                auto result = sub(n1, shallow_copy);
                shallow_copy.byte = nullptr;
                Q_UNUSED(shallow_copy)
                return result;
            }
            else {
                Scalar shallow_copy(const_cast<ScalarUnit*>(n1.byte), -n1.length, n1.power);
                auto result = sub(n2, shallow_copy);
                shallow_copy.byte = nullptr;
                Q_UNUSED(shallow_copy)
                return result;
            }
        }
        else {
            const Scalar* big;
            const Scalar* small;
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
            //Estimate the ed of result first, will calculate it accurately later.
            int length = (big->power + std::max(bigSize - big->power, smallSize - small->power));
            length = length > BasicConst::getInstance().GlobalPrecision
                     ? BasicConst::getInstance().GlobalPrecision : length;
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
                //useableSmall is the part whose add result will fall in GlobalPrecision.
                int useableSmall = small->power - (big->power - BasicConst::getInstance().GlobalPrecision);
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
    inline Scalar Scalar::sub(const Scalar& n1, const Scalar& n2) {
        if(n1.isZero())
            return -n2;
        else if(n2.isZero())
            return Scalar(n1);
        else if (n1.length > 0) {
            if (n2.length < 0) {
                Scalar shallow_copy(const_cast<ScalarUnit*>(n2.byte), -n2.length, n2.power);
                Scalar result = add(n1, shallow_copy);
                shallow_copy.byte = nullptr;
                Q_UNUSED(shallow_copy)
                return result;
            }
            else {
                const Scalar* big;
                const Scalar* small;
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
                //Estimate the ed of result first, will calculate it accurately later.
                int length = (big->power + std::max(bigSize - big->power, smallSize - small->power));
                length = length > BasicConst::getInstance().GlobalPrecision
                         ? BasicConst::getInstance().GlobalPrecision : length;
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
                    //useableSmall is the part whose sub result will fall in GlobalPrecision.
                    int useableSmall = small->power - (big->power - BasicConst::getInstance().GlobalPrecision);
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
            Scalar shallow_copy(const_cast<ScalarUnit*>(n1.byte), -n1.length, n1.power);
            if (n2.length > 0) {
                Scalar result = add(shallow_copy, n2);
                result.toOpposite();
                shallow_copy.byte = nullptr;
                Q_UNUSED(shallow_copy)
                return result;
            }
            else {
                Scalar shallow_copy_1(const_cast<ScalarUnit*>(n2.byte), n2.length, n2.power);
                Scalar result = sub(shallow_copy_1, shallow_copy);
                shallow_copy.byte = shallow_copy_1.byte = nullptr;
                Q_UNUSED(shallow_copy)
                Q_UNUSED(shallow_copy_1)
                return result;
            }
        }
    }

    inline Scalar Scalar::mul(const Scalar& n1, const Scalar& n2) {
        if(n1 == BasicConst::getInstance().get_1())
            return Scalar(n2);
        else if(n2 == BasicConst::getInstance().get_1())
            return Scalar(n1);
        else {
            const int size1 = n1.getSize();
            const int size2 = n2.getSize();
            //Estimate the ed of result first. we will calculate it accurately later.
            auto length = size1 + size2;
            auto byte = reinterpret_cast<ScalarUnit*>(calloc(length, sizeof(ScalarUnit)));
            for (int i = 0; i < size1; ++i)
                byte[i + size2] = mulAddArrByWord(byte + i, n2.byte, size2, n1.byte[i]);
            ///////////////////////////////////////Get byte, length and power//////////////////////////;
            int power = n1.power + n2.power + 1;
            if (byte[length - 1] == 0) {
                --length;
                --power;
                byte = reinterpret_cast<ScalarUnit*>(realloc(byte, length * sizeof(ScalarUnit)));
            }
            ////////////////////////////////////Out put////////////////////////////////////////
            return Scalar(byte, (n1.length ^ n2.length) < 0 ? -length : length, power); //NOLINT
        }
    }

    inline Scalar Scalar::div(const Scalar& n1, const Scalar& n2) {
        if(Q_UNLIKELY(n2.isZero()))
            qFatal("Encountered dividing by zero exception.");

        if(!n1.isZero()) {
            if(n2 != BasicConst::getInstance().get_1()) {
                const auto n1_size = n1.getSize(), n2_size = n2.getSize();
                //Add one to avoid precision loss during right shift.
                auto arr1_len = std::max(n1_size, n2_size) + 1;
                auto n1_blank = arr1_len - n1_size;
                auto arr1 = new ScalarUnit[arr1_len];
                memcpy(arr1 + n1_blank, n1.byte, n1_size * sizeof(ScalarUnit));
                memset(arr1, 0, n1_blank * sizeof(ScalarUnit));
                //Size of arr2 is arranged 1 less than arr1.
                auto arr2_len = arr1_len - 1;
                auto n2_blank = arr2_len - n2_size;
                auto arr2 = new ScalarUnit[arr2_len];
                memcpy(arr2 + n2_blank, n2.byte, n2_size * sizeof(ScalarUnit));
                memset(arr2, 0, n2_blank * sizeof(ScalarUnit));

                const int n1_shift = static_cast<int>(countLeadingZeros(n1.byte[n1_size - 1])) - 1;
                if(n1_shift > 0)
                    byteLeftShiftEq(arr1, arr1_len, n1_shift);
                else
                    byteRightShiftEq(arr1, arr1_len, -n1_shift);
                const int n2_shift = static_cast<int>(countLeadingZeros(n2.byte[n2_size - 1]));
                byteLeftShiftEq(arr2, arr2_len, n2_shift);
                ////////////////////////////////Calculate cursory first//////////////////////////////////////
                //Estimate the length of result.
                int length = BasicConst::getInstance().GlobalPrecision;
                //let n1_copy's power equal to n2_copy, power of the result will change correspondingly.
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
                return Scalar(byte, (n1.length ^ n2.length) < 0 ? -length : length //NOLINT
                        , n1.getPower() - n2.getPower() - 1, 1) >> (n1_shift - n2_shift);
            }
            else
                return Scalar(n1);
        }
        else
            return getZero();
    }
    /*
     * If the length of new array is larger than GlobalPrecision, it will be set to GlobalPrecision.
     * Return true if array is cut.
     */
    inline bool Scalar::cutLength(Scalar& n) {
        bool result = false;
        int size = n.getSize();

        if(size > BasicConst::getInstance().GlobalPrecision) {
            result = true;
            int cutFrom = size - BasicConst::getInstance().GlobalPrecision;
            auto new_byte = reinterpret_cast<ScalarUnit*>(malloc(BasicConst::getInstance().GlobalPrecision * sizeof(ScalarUnit)));
            memcpy(new_byte, n.byte + cutFrom, BasicConst::getInstance().GlobalPrecision * sizeof(ScalarUnit));
            free(n.byte);
            n.byte = new_byte;
            auto length = BasicConst::getInstance().GlobalPrecision;
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
    inline void Scalar::cutZero(Scalar& n) {
        int size = n.getSize();
        int id = size - 1;
        while(n.byte[id] == 0 && id > 0)
            --id;
        ++id;

        if(id != size) {
            int shorten = size - id;
            n.byte = reinterpret_cast<ScalarUnit*>(realloc(n.byte, id * sizeof(ScalarUnit)));
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
}

#endif
