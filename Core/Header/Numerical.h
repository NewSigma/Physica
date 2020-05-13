/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_NUMERICAL_H
#define PHYSICA_NUMERICAL_H

#include <iostream>
#include "Const.h"
#include "ElementaryFunction.h"
#include "SystemBits.h"

namespace Physica::Core {
    class Numerical {
        //Store effective digits.
        NumericalUnit* byte;
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
        Numerical(NumericalUnit*& byte, int length, int power, NumericalUnit a = 0) noexcept;
        Numerical(NumericalUnit*&& byte, int length, int power, NumericalUnit a = 0) noexcept;
        Numerical(const Numerical& n) noexcept;
        Numerical(Numerical&& n) noexcept;
        explicit Numerical(const Numerical* n) noexcept;
        explicit Numerical(SignedNumericalUnit unit, NumericalUnit a = 0) noexcept;
        explicit Numerical(double d, NumericalUnit a = 0);
        explicit Numerical(const char* s, NumericalUnit a = 0);
        explicit Numerical(const wchar_t* s, NumericalUnit a = 0);
        explicit Numerical(const std::string& s, NumericalUnit a = 0);
        explicit Numerical(const std::wstring& s, NumericalUnit a = 0);
        ~Numerical();

        explicit operator double() const;
        friend std::ostream& operator<<(std::ostream& os, const Numerical& n);
        Numerical operator<<(int bits);
        Numerical operator>>(int bits);
        inline void operator<<=(int bits) noexcept { power += bits; }
        inline void operator>>=(int bits) noexcept { power -= bits; }
        inline NumericalUnit& operator[](unsigned int index) { return byte[index]; }
        inline const NumericalUnit& operator[](unsigned int index) const { return byte[index]; }
        Numerical& operator=(const Numerical& n);
        Numerical& operator=(Numerical&& n) noexcept;
        Numerical operator^(const Numerical& n) const;
        Numerical operator-() const;

        Numerical& applyError(const Numerical& error);

        [[nodiscard]] int getLength() const noexcept { return length; }
        [[nodiscard]] int getPower() const noexcept { return power; }
        [[nodiscard]] NumericalUnit getA() const noexcept { return a; }
        [[nodiscard]] int getSize() const noexcept { return abs(length); }
        Numerical& toAbs() noexcept { length = getSize(); return *this; } //NOLINT
        Numerical& toOpposite() noexcept { length = -length; return *this; }
        Numerical& toUnitA() noexcept { a = 1; return *this; }
        Numerical& clearA() noexcept { a = 0; return *this; }

        [[nodiscard]] bool isZero() const { return byte[getSize() - 1] == 0; }
        [[nodiscard]] bool isPositive() const { return !isZero() && length > 0; }
        [[nodiscard]] bool isNegative() const { return !isZero() && length < 0; }
        [[nodiscard]] bool isInteger() const { return getSize() == power + 1; }

        friend Numerical add (const Numerical& n1, const Numerical& n2);
        friend Numerical sub (const Numerical& n1, const Numerical& n2);
        friend Numerical mul (const Numerical& n1, const Numerical& n2);
        friend Numerical div (const Numerical& n1, const Numerical& n2);
        friend Numerical square(const Numerical& n);
        friend bool cutLength(Numerical& n);
        friend void cutZero(Numerical& n);
        friend Numerical sqrt_light(const Numerical& n);
    };
////////////////////////////////Operators/////////////////////////////////////
    Numerical operator+(const Numerical& n1, const Numerical& n2);
    Numerical operator-(const Numerical& n1, const Numerical& n2);
    Numerical operator*(const Numerical& n1, const Numerical& n2);
    Numerical operator/(const Numerical& n1, const Numerical& n2);
    inline void operator+=(Numerical& n1, const Numerical& n2) { n1 = n1 + n2; }
    inline void operator-=(Numerical& n1, const Numerical& n2) { n1 = n1 - n2; }
    inline void operator*=(Numerical& n1, const Numerical& n2) { n1 = n1 * n2; }
    inline void operator/=(Numerical& n1, const Numerical& n2) { n1 = n1 / n2; }
    inline void operator^=(Numerical& n1, const Numerical& n2) { n1 = n1 ^ n2; }
    bool operator>(const Numerical& n1, const Numerical& n2);
    bool operator<(const Numerical& n1, const Numerical& n2);
    bool operator>=(const Numerical& n1, const Numerical& n2);
    bool operator<=(const Numerical& n1, const Numerical& n2);
    bool operator==(const Numerical& n1, const Numerical& n2);
    bool operator!=(const Numerical& n1, const Numerical& n2);
    inline Numerical& operator++(Numerical& n) { n += BasicConst::getInstance().get_1(); return n; }
    inline Numerical& operator--(Numerical& n) { n -= BasicConst::getInstance().get_1(); return n; }
    inline Numerical operator++(Numerical& n, int) { Numerical temp(n); n += BasicConst::getInstance().get_1(); return temp; } //NOLINT
    inline Numerical operator--(Numerical& n, int) { Numerical temp(n); n -= BasicConst::getInstance().get_1(); return temp; } //NOLINT
////////////////////////////////Helper functions/////////////////////////////////////
    Numerical getAccuracy(const Numerical& n);
    inline Numerical getMaximum(const Numerical& n) { return add(n, getAccuracy(n)).clearA(); }
    inline Numerical getMinimum(const Numerical& n) { return sub(n, getAccuracy(n)).clearA(); }
    inline Numerical getZero() { return Numerical(BasicConst::getInstance().get_0()); }
    inline Numerical getOne() { return Numerical(BasicConst::getInstance().get_1()); }
    inline Numerical getTwo() { return Numerical(BasicConst::getInstance().get_2()); }
}

#endif
