/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_NUMERICAL_H
#define PHYSICA_NUMERICAL_H

#include <iostream>
#include "Const.h"
#include "CalcBasic.h"
#include "ElementaryFunction.h"

extern const BasicConst* basicConst;

class Numerical {
    //Store effective digits.
    //TODO Depending on the platform, unsigned long may not be the best choice but the type whose length equals to the word length is.
    unsigned long* byte;
    //Number = (x0 +- a * 10^(1-length)) * 10^power
    int power;
    /*
     * Length of byte = abs(length).
     * sign of length and sign of Numerical are same. (when Numerical != 0)
     */
    int length;
    //Accuracy
    unsigned long a;
public:
    Numerical(unsigned long* byte, int length, int power, unsigned long a = 0);
    Numerical(const Numerical& n);
    explicit Numerical(const Numerical* n);
    explicit Numerical(double d, unsigned long a = 0);
    explicit Numerical(const char* s, unsigned long a = 0);
    explicit Numerical(const wchar_t* s, unsigned long a = 0);
    explicit Numerical(const std::string& s, unsigned long a = 0);
    explicit Numerical(const std::wstring& s, unsigned long a = 0);
    ~Numerical();

    explicit operator double() const;
    friend std::ostream& operator<<(std::ostream& os, const Numerical& n);
    void operator<<(int bits);
    void operator>>(int bits);
    unsigned long operator[](unsigned int index);
    Numerical& operator= (const Numerical& n);
    Numerical* operator+ (const Numerical& n) const;
    Numerical* operator- (const Numerical& n) const;
    Numerical* operator* (const Numerical& n) const;
    Numerical* operator/ (const Numerical& n) const;
    Numerical* operator^ (const Numerical& n) const;
    void operator+= (const Numerical& n);
    void operator-= (const Numerical& n);
    void operator*= (const Numerical& n);
    void operator/= (const Numerical& n);
    void operator^= (const Numerical& n);
    bool operator> (const Numerical& n) const;
    bool operator< (const Numerical& n) const;
    bool operator>= (const Numerical& n) const;
    bool operator<= (const Numerical& n) const;
    bool operator== (const Numerical& n) const;
    bool operator!= (const Numerical& n) const;
    Numerical* operator- () const;

    Numerical* getAccuracy() const;
    Numerical* getMaximum() const;
    Numerical* getMinimum() const;
    void applyError(const Numerical* error);

    const int& getLength() const { return length; }
    const int& getPower() const { return power; }
    const unsigned long& getA() const { return a; }
    int getSize() const { return abs(length); }
    void toAbs() { length = getSize(); }
    void toOpposite() { length = -length; }
    void toUnitA() { a = 1; }
    void clearA() { a = 0; }

    bool isZero() const { return byte[getSize() - 1] == 0; }
    bool isPositive() const { return !isZero() && length > 0; }
    bool isNegative() const { return !isZero() && length < 0; }
    bool isInteger() const { return getSize() == power + 1; }

    friend Numerical* add (const Numerical& n1, const Numerical& n2);
    friend Numerical* subtract (const Numerical& n1, const Numerical& n2);
    friend Numerical* multiply (const Numerical& n1, const Numerical& n2);
    friend Numerical* divide (const Numerical& n1, const Numerical& n2);
    friend bool cutLength(Numerical* n);
    friend void cutZero(Numerical* n);
    friend Numerical* sqrt_light(const Numerical& n);
};
////////////////////////////////Helper functions/////////////////////////////////////
inline void move(Numerical*& dest, Numerical*&& src) { delete dest; dest = src; }
inline Numerical* getZero() { return new Numerical(basicConst->get_0()); }
inline Numerical* getOne() { return new Numerical(basicConst->get_1()); }
inline Numerical* getTwo() { return new Numerical(basicConst->get_2()); }

#endif
