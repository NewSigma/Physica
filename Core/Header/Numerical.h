/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_NUMERICAL_H
#define PHYSICA_NUMERICAL_H

#include <iostream>
#include "Const.h"
#include "RealNum.h"

extern const Const_1* const_1;

class Numerical {
public:
    //Store effective digits.
    unsigned char* byte;
    //Number = (x0 +- a * 10^(1-length)) * 10^power
    int power;
    /*
     * Length of byte = abs(length).
     * sign of length and sign of Numerical are same. (when Numerical != 0)
     */
    signed char length;
    //Accuracy
    unsigned char a = 0;

    Numerical();
    Numerical(const Numerical& n);
    Numerical(unsigned char* byte, signed char length, int power, unsigned char acc = 0);
    explicit Numerical(const Numerical* n);
    explicit Numerical(double d, unsigned char acc = 0);
    explicit Numerical(const char* s, unsigned char acc = 0);
    explicit Numerical(const wchar_t* s, unsigned char acc = 0);
    explicit Numerical(const std::string& s, unsigned char acc = 0);
    explicit Numerical(const std::wstring& s, unsigned char acc = 0);
    ~Numerical();

    std::string toString() const;
    explicit operator double() const;

    friend std::ostream& operator<<(std::ostream& os, const Numerical& n);
    void operator<<(Numerical& n);
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
    bool applyError(const Numerical* copy_error);
    inline int getSize() const { return abs(length); }
    inline bool isZero() const { return byte[getSize() - 1] == 0; }
    inline bool isPositive() const { return !isZero() && length > 0; }
    inline bool isNegative() const { return !isZero() && length < 0; }
    inline bool isInteger() const { return length == power + 1; }
};
////////////////////////////////Helper functions/////////////////////////////////////
void reverse(unsigned char* arr, int length);
inline Numerical* getZero() { return new Numerical(const_1->_0); }
inline Numerical* getOne() { return new Numerical(const_1->_1); }
inline Numerical* getTwo() { return new Numerical(const_1->_2); }
//////////////////////////////Process functions////////////////////////////////////////
Numerical* add (const Numerical& n1, const Numerical& n2);
Numerical* subtract (const Numerical& n1, const Numerical& n2);
Numerical* multiply (const Numerical& n1, const Numerical& n2);
Numerical* divide (const Numerical& n1, const Numerical& n2);
bool cutLength(Numerical* n);
void cutZero(Numerical* n);
////////////////////////////////////////Elementary Functions////////////////////////////////////////////
Numerical* randomNumerical();
Numerical* randomNumerical(Numerical* lowerBound, Numerical* upperBound);
Numerical* reciprocal(const Numerical& n);
Numerical* sqrt_light(const Numerical& n);
Numerical* sqrt(const Numerical& n);
Numerical* factorial(const Numerical& n);
Numerical* ln_light(const Numerical& n);
Numerical* ln(const Numerical& n);
Numerical* log(const Numerical& n, const Numerical& a);
Numerical* exp(const Numerical& n);
Numerical* pow(const Numerical& n, const Numerical& a);
Numerical* cos(const Numerical& n);
Numerical* sin(const Numerical& n);
Numerical* tan(const Numerical& n);
Numerical* sec(const Numerical& n);
Numerical* csc(const Numerical& n);
Numerical* cot(const Numerical& n);
Numerical* arccos(const Numerical& n);
Numerical* arcsin(const Numerical& n);
Numerical* arctan(const Numerical& n);
Numerical* arcsec(const Numerical& n);
Numerical* arccsc(const Numerical& n);
Numerical* arccot(const Numerical& n);
Numerical* cosh(const Numerical& n);
Numerical* sinh(const Numerical& n);
Numerical* tanh(const Numerical& n);
Numerical* sech(const Numerical& n);
Numerical* csch(const Numerical& n);
Numerical* coth(const Numerical& n);
Numerical* arccosh(const Numerical& n);
Numerical* arcsinh(const Numerical& n);
Numerical* arctanh(const Numerical& n);
Numerical* arcsech(const Numerical& n);
Numerical* arccsch(const Numerical& n);
Numerical* arccoth(const Numerical& n);

#endif
