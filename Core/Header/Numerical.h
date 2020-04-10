/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_NUMERICAL_H
#define PHYSICA_NUMERICAL_H

#include <iostream>
#include <climits>
#include "Const.h"
#include "RealNum.h"

extern const BasicConst* basicConst;

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
    unsigned char a;

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
    inline bool isInteger() const { return getSize() == power + 1; }
};
////////////////////////////////Helper functions/////////////////////////////////////
inline Numerical* getZero() { return new Numerical(basicConst->get_0()); }
inline Numerical* getOne() { return new Numerical(basicConst->get_1()); }
inline Numerical* getTwo() { return new Numerical(basicConst->get_2()); }

#endif
