/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_NUMERICAL_H
#define PHYSICA_NUMERICAL_H

#include <iostream>
#include <cstring>
#include <QtCore/qlogging.h>
#include "Core/Header/SystemBits.h"
#include "Core/Header/Const.h"
#include "Core/MultiplePrecision/Basic/DivBasic.h"
#include "Core/MultiplePrecision/Util/ArraySupport.h"
#include "Core/Header/ElementaryFunction.h"

namespace Physica::Core {
    class Numerical {
        //Store effective digits.
        NumericalUnit* __restrict byte;
        /*
         * Length of byte = abs(length).
         * sign of length and sign of Numerical are same. (when Numerical != 0)
        */
        int length;
        /*
         * Number = (x0 +- a * (2 ^ __WORDSIZE) ^ (1-length)) * (2 ^ __WORDSIZE) ^power
         * We have not considered overflow of power in our codes elsewhere.
         */
        int power;
        //Accuracy
        NumericalUnit a;
    public:
        Numerical() noexcept;
        Numerical(int length, int power, NumericalUnit a = 0) noexcept;
        Numerical(const Numerical& n) noexcept;
        Numerical(Numerical&& n) noexcept;
        explicit Numerical(SignedNumericalUnit unit, NumericalUnit a = 0) noexcept;
        explicit Numerical(double d, NumericalUnit a = 0);
        explicit Numerical(const char* s, NumericalUnit a = 0);
        explicit Numerical(const wchar_t* s, NumericalUnit a = 0);
        explicit Numerical(const std::string& s, NumericalUnit a = 0);
        explicit Numerical(const std::wstring& s, NumericalUnit a = 0);
        ~Numerical();
        /* Operators */
        explicit operator double() const;
        friend std::ostream& operator<<(std::ostream& os, const Numerical& n);
        Numerical operator+(const Numerical& n) const;
        Numerical operator-(const Numerical& n) const;
        Numerical operator*(const Numerical& n) const;
        Numerical operator/(const Numerical& n) const;
        Numerical operator<<(int bits) const;
        Numerical operator>>(int bits) const;
        void operator<<=(int bits) noexcept { *this = *this << bits; }
        void operator>>=(int bits) noexcept { *this = *this >> bits; }
        NumericalUnit& operator[](unsigned int index) { return byte[index]; }
        const NumericalUnit& operator[](unsigned int index) const { return byte[index]; }
        Numerical& operator=(const Numerical& n);
        Numerical& operator=(Numerical&& n) noexcept;
        Numerical operator^(const Numerical& n) const;
        Numerical operator-() const;
        /* Helpers */
        Numerical& applyError(const Numerical& error);
        void swap(Numerical& n) noexcept;
        Numerical& toAbs() noexcept { length = getSize(); return *this; }
        Numerical& toOpposite() noexcept { length = -length; return *this; }
        Numerical& toUnitA() noexcept { a = 1; return *this; }
        Numerical& clearA() noexcept { a = 0; return *this; }
        /* Getters */
        [[nodiscard]] int getLength() const noexcept { return length; }
        [[nodiscard]] int getPower() const noexcept { return power; }
        [[nodiscard]] NumericalUnit getA() const noexcept { return a; }
        [[nodiscard]] int getSize() const noexcept { return abs(length); }
        [[nodiscard]] bool isZero() const { return byte[getSize() - 1] == 0; }
        [[nodiscard]] bool isPositive() const { return !isZero() && length > 0; }
        [[nodiscard]] bool isNegative() const { return !isZero() && length < 0; }
        [[nodiscard]] bool isInteger() const { return getSize() == power + 1; }
        [[nodiscard]] Numerical getAccuracy() const;
        [[nodiscard]] Numerical getMaximum() const { return Numerical::add(*this, getAccuracy()).clearA(); }
        [[nodiscard]] Numerical getMinimum() const { return Numerical::sub(*this, getAccuracy()).clearA(); }
    protected:
        Numerical(NumericalUnit* byte, int length, int power, NumericalUnit a = 0);
        inline static Numerical add (const Numerical& n1, const Numerical& n2);
        inline static Numerical sub (const Numerical& n1, const Numerical& n2);
        inline static Numerical mul (const Numerical& n1, const Numerical& n2);
        inline static Numerical div (const Numerical& n1, const Numerical& n2);
        inline static bool cutLength(Numerical& n);
        inline static void cutZero(Numerical& n);
        /* Friends */
        friend class Solve;
        friend Numerical square(const Numerical& n);
        friend Numerical sqrt_light(const Numerical& n);
        friend Numerical ln_light(const Numerical& n);
    };
    /* Operators */
    bool operator>(const Numerical& n1, const Numerical& n2);
    bool operator<(const Numerical& n1, const Numerical& n2);
    bool operator>=(const Numerical& n1, const Numerical& n2);
    bool operator<=(const Numerical& n1, const Numerical& n2);
    bool operator==(const Numerical& n1, const Numerical& n2);
    bool operator!=(const Numerical& n1, const Numerical& n2);
    /* Inline Implementations */
    inline void operator+=(Numerical& n1, const Numerical& n2) { n1 = n1 + n2; }
    inline void operator-=(Numerical& n1, const Numerical& n2) { n1 = n1 - n2; }
    inline void operator*=(Numerical& n1, const Numerical& n2) { n1 = n1 * n2; }
    inline void operator/=(Numerical& n1, const Numerical& n2) { n1 = n1 / n2; }
    inline void operator^=(Numerical& n1, const Numerical& n2) { n1 = n1 ^ n2; }
    inline Numerical& operator++(Numerical& n) { n += BasicConst::getInstance().get_1(); return n; }
    inline Numerical& operator--(Numerical& n) { n -= BasicConst::getInstance().get_1(); return n; }
    inline Numerical operator++(Numerical& n, int) { Numerical temp(n); n += BasicConst::getInstance().get_1(); return temp; } //NOLINT
    inline Numerical operator--(Numerical& n, int) { Numerical temp(n); n -= BasicConst::getInstance().get_1(); return temp; } //NOLINT
    inline Numerical getZero() { return Numerical(BasicConst::getInstance().get_0()); }
    inline Numerical getOne() { return Numerical(BasicConst::getInstance().get_1()); }
    inline Numerical getTwo() { return Numerical(BasicConst::getInstance().get_2()); }
    inline void swap(Numerical& n1, Numerical& n2) noexcept { n1.swap(n2); }
    /*
     * The following four functions simply calculate the result while operator functions will
     * consider the accuracy.
     */
    //Refactor: move the key algorithm into a separated method. We may not use every words.
    inline Numerical Numerical::add(const Numerical& n1, const Numerical& n2) {
        if(n1.isZero())
            return Numerical(n2);
        else if(n2.isZero())
            return Numerical(n1);
        else if ((n1.length ^ n2.length) < 0) { // NOLINT(hicpp-signed-bitwise)
            if (n1.length > 0) {
                Numerical shallow_copy(const_cast<NumericalUnit*>(n2.byte), -n2.length, n2.power);
                auto result = sub(n1, shallow_copy);
                shallow_copy.byte = nullptr;
                Q_UNUSED(shallow_copy)
                return result;
            }
            else {
                Numerical shallow_copy(const_cast<NumericalUnit*>(n1.byte), -n1.length, n1.power);
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
            auto byte = reinterpret_cast<NumericalUnit*>(malloc(length * sizeof(NumericalUnit)));
            memcpy(byte + length - bigSize, big->byte, bigSize * sizeof(NumericalUnit));
            memset(byte, 0, (length - bigSize) * sizeof(NumericalUnit));

            int index = length - big->power + small->power - smallSize;
            NumericalUnit aByte;
            NumericalUnit carry = 0;
            NumericalUnit carry_temp;
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
                byte = reinterpret_cast<NumericalUnit*>(realloc(byte, length * sizeof(NumericalUnit)));
                byte[length - 1] = 1;
            }
            ////////////////////////////////////Out put////////////////////////////////////////
            return Numerical(byte, big->length < 0 ? -length : length, power);
        }
    }
    //Refactor: move the key algorithm into a separated method. We may not use every words.
    inline Numerical Numerical::sub(const Numerical& n1, const Numerical& n2) {
        if(n1.isZero())
            return -n2;
        else if(n2.isZero())
            return Numerical(n1);
        else if (n1.length > 0) {
            if (n2.length < 0) {
                Numerical shallow_copy(const_cast<NumericalUnit*>(n2.byte), -n2.length, n2.power);
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
                auto byte = reinterpret_cast<NumericalUnit*>(malloc(length * sizeof(NumericalUnit)));
                memcpy(byte + length - bigSize, big->byte, bigSize * sizeof(NumericalUnit));
                memset(byte, 0, (length - bigSize) * sizeof(NumericalUnit));

                int index = length - big->power + small->power - smallSize;
                NumericalUnit aByte;
                NumericalUnit carry = 0;
                NumericalUnit carry_temp;
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
                    free(byte);
                    goto redo;
                }
                Numerical result(byte, changeSign ? -length : length, big->power);
                cutZero(result);
                return result;
            }
        }
        else {
            Numerical shallow_copy(const_cast<NumericalUnit*>(n1.byte), -n1.length, n1.power);
            if (n2.length > 0) {
                Numerical result = add(shallow_copy, n2);
                result.toOpposite();
                shallow_copy.byte = nullptr;
                Q_UNUSED(shallow_copy)
                return result;
            }
            else {
                Numerical shallow_copy_1(const_cast<NumericalUnit*>(n2.byte), n2.length, n2.power);
                Numerical result = sub(shallow_copy_1, shallow_copy);
                shallow_copy.byte = shallow_copy_1.byte = nullptr;
                Q_UNUSED(shallow_copy)
                Q_UNUSED(shallow_copy_1)
                return result;
            }
        }
    }

    inline Numerical Numerical::mul(const Numerical& n1, const Numerical& n2) {
        if(n1 == BasicConst::getInstance().get_1())
            return Numerical(n2);
        else if(n2 == BasicConst::getInstance().get_1())
            return Numerical(n1);
        else {
            const int size1 = n1.getSize();
            const int size2 = n2.getSize();
            //Estimate the ed of result first. we will calculate it accurately later.
            auto length = size1 + size2;
            auto byte = reinterpret_cast<NumericalUnit*>(calloc(length, sizeof(NumericalUnit)));
            for (int i = 0; i < size1; ++i)
                byte[i + size2] = mulAddArrByWord(byte + i, n2.byte, size2, n1.byte[i]);
            ///////////////////////////////////////Get byte, length and power//////////////////////////;
            int power = n1.power + n2.power + 1;
            if (byte[length - 1] == 0) {
                --length;
                --power;
                byte = reinterpret_cast<NumericalUnit*>(realloc(byte, length * sizeof(NumericalUnit)));
            }
            ////////////////////////////////////Out put////////////////////////////////////////
            return Numerical(byte, (n1.length ^ n2.length) < 0 ? -length : length, power); //NOLINT
        }
    }

    inline Numerical Numerical::div(const Numerical& n1, const Numerical& n2) {
        if(Q_UNLIKELY(n2.isZero()))
            qFatal("Encountered dividing by zero exception.");

        if(!n1.isZero()) {
            if(n2 != BasicConst::getInstance().get_1()) {
                auto n1_size = n1.getSize(), n2_size = n2.getSize();
                //Add one to avoid precision loss during right shift.
                auto arr1_len = std::max(n1_size, n2_size) + 1;
                auto n1_blank = arr1_len - n1_size;
                auto arr1 = new NumericalUnit[arr1_len];
                memcpy(arr1 + n1_blank, n1.byte, n1_size * sizeof(NumericalUnit));
                memset(arr1, 0, n1_blank * sizeof(NumericalUnit));
                //Size of arr2 is arranged 1 less than arr1.
                auto arr2_len = arr1_len - 1;
                auto n2_blank = arr2_len - n2_size;
                auto arr2 = new NumericalUnit[arr2_len];
                memcpy(arr2 + n2_blank, n2.byte, n2_size * sizeof(NumericalUnit));
                memset(arr2, 0, n2_blank * sizeof(NumericalUnit));

                int n1_shift = static_cast<int>(countLeadingZeros(n1.byte[n1.getSize() - 1])) - 1;
                if(n1_shift > 0)
                    byteLeftShiftEq(arr1, arr1_len, n1_shift);
                else
                    byteRightShiftEq(arr1, arr1_len, -n1_shift);
                int n2_shift = static_cast<int>(countLeadingZeros(n2.byte[n2.getSize() - 1]));
                byteLeftShiftEq(arr2, arr2_len, n2_shift);
                ////////////////////////////////Calculate cursory first//////////////////////////////////////
                //Estimate the length of result.
                int length = BasicConst::getInstance().GlobalPrecision;
                //let n1_copy's power equal to n2_copy, power of the result will change correspondingly.
                auto byte = reinterpret_cast<NumericalUnit*>(malloc(length * sizeof(NumericalUnit)));

                for(int i = length - 1; i >= 0; --i) {
                    byte[i] = divArrByFullArrWith1Word(arr1, arr2, arr2_len);
                    arr1[arr2_len] -= mulSubArrByWord(arr1, arr2, arr2_len, byte[i]);
                    byteLeftShiftEq(arr1, arr1_len, NumericalUnitWidth);
                }
                delete[] arr1;
                delete[] arr2;
                ////////////////////////////////////Out put////////////////////////////////////////
                //1 comes from the algorithm
                return Numerical(byte, (n1.length ^ n2.length) < 0 ? -length : length //NOLINT
                        , n1.getPower() - n2.getPower() - 1, 1) >> (n1_shift - n2_shift);
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
    inline bool Numerical::cutLength(Numerical& n) {
        bool result = false;
        int size = n.getSize();

        if(size > BasicConst::getInstance().GlobalPrecision) {
            result = true;
            int cutFrom = size - BasicConst::getInstance().GlobalPrecision;
            auto new_byte = reinterpret_cast<NumericalUnit*>(malloc(BasicConst::getInstance().GlobalPrecision * sizeof(NumericalUnit)));
            memcpy(new_byte, n.byte + cutFrom, BasicConst::getInstance().GlobalPrecision * sizeof(NumericalUnit));
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
    inline void Numerical::cutZero(Numerical& n) {
        int size = n.getSize();
        int id = size - 1;
        while(n.byte[id] == 0 && id > 0)
            --id;
        ++id;

        if(id != size) {
            int shorten = size - id;
            n.byte = reinterpret_cast<NumericalUnit*>(realloc(n.byte, id * sizeof(NumericalUnit)));
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
