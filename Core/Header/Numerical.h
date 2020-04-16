/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_NUMERICAL_H
#define PHYSICA_NUMERICAL_H

#include <iostream>
#include "Const.h"
#include "CalcBasic.h"
#include "ElementaryFunction.h"
#include "SystemBits.h"

extern const BasicConst* basicConst;

class Numerical {
    //Store effective digits.
    NumericalUnit* byte;
    /*
     * Length of byte = abs(length).
     * sign of length and sign of Numerical are same. (when Numerical != 0)
    */
    int length;
    //Number = (x0 +- a * 10^(1-length)) * 10^power
    int power;
    //Accuracy
    NumericalUnit a;
public:
    Numerical(NumericalUnit* byte, int length, int power, NumericalUnit a = 0);
    Numerical(const Numerical& n);
    Numerical(Numerical&& n) noexcept;
    explicit Numerical(const Numerical* n);
    explicit Numerical(double d, NumericalUnit a = 0);
    explicit Numerical(const char* s, NumericalUnit a = 0);
    explicit Numerical(const wchar_t* s, NumericalUnit a = 0);
    explicit Numerical(const std::string& s, NumericalUnit a = 0);
    explicit Numerical(const std::wstring& s, NumericalUnit a = 0);
    ~Numerical();

    explicit operator double() const;
    friend std::ostream& operator<<(std::ostream& os, const Numerical& n);
    void operator<<(int bits);
    void operator>>(int bits);
    NumericalUnit operator[](unsigned int index) const;
    Numerical& operator=(const Numerical& n);
    Numerical& operator=(Numerical&& n) noexcept;
    Numerical operator+(const Numerical& n) const;
    Numerical operator-(const Numerical& n) const;
    Numerical operator*(const Numerical& n) const;
    Numerical operator/(const Numerical& n) const;
    Numerical operator^(const Numerical& n) const;
    void operator+=(const Numerical& n) { *this = *this + n; };
    void operator-=(const Numerical& n) { *this = *this - n; };
    void operator*=(const Numerical& n) { *this = *this * n; };
    void operator/=(const Numerical& n) { *this = *this / n; };
    void operator^=(const Numerical& n) { *this = *this ^ n; };
    bool operator>(const Numerical& n) const;
    bool operator<(const Numerical& n) const;
    bool operator>=(const Numerical& n) const;
    bool operator<=(const Numerical& n) const;
    bool operator==(const Numerical& n) const;
    bool operator!=(const Numerical& n) const;
    Numerical operator-() const;

    Numerical getAccuracy() const;
    Numerical getMaximum() const { return add(*this, getAccuracy()); }
    Numerical getMinimum() const { return sub(*this, getAccuracy()); }
    Numerical& applyError(const Numerical& error);

    const int& getLength() const { return length; }
    const int& getPower() const { return power; }
    const NumericalUnit& getA() const { return a; }
    int getSize() const { return abs(length); }
    Numerical& toAbs() { length = getSize(); return *this; }
    Numerical& toOpposite() { length = -length; return *this; }
    Numerical& toUnitA() { a = 1; return *this; }
    Numerical& clearA() { a = 0; return *this; }

    bool isZero() const { return byte[getSize() - 1] == 0; }
    bool isPositive() const { return !isZero() && length > 0; }
    bool isNegative() const { return !isZero() && length < 0; }
    bool isInteger() const { return getSize() == power + 1; }

    friend Numerical add (const Numerical& n1, const Numerical& n2);
    friend Numerical sub (const Numerical& n1, const Numerical& n2);
    friend Numerical mul (const Numerical& n1, const Numerical& n2);
    friend Numerical div (const Numerical& n1, const Numerical& n2);
    friend bool cutLength(Numerical& n);
    friend void cutZero(Numerical& n);
    friend Numerical sqrt_light(const Numerical& n);
};
////////////////////////////////Helper functions/////////////////////////////////////
void printElements(const Numerical& n);
inline Numerical getZero() { return Numerical(basicConst->get_0()); }
inline Numerical getOne() { return Numerical(basicConst->get_1()); }
inline Numerical getTwo() { return Numerical(basicConst->get_2()); }

#endif
