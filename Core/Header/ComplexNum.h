/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_COMPLEXNUM_H
#define PHYSICA_COMPLEXNUM_H

#include "AbstractNum.h"
class Numerical;
class Vector;

class ComplexNum : public AbstractNum{
public:
    Numerical* real;
    Numerical* imagine;

    ComplexNum(Numerical* n1, Numerical* n2, bool polar = false);
    ComplexNum(ComplexNum& instance);
    ComplexNum(ComplexNum* instance);
    virtual ~ComplexNum();
    ComplexNum* toConjugate() const;
    Numerical* toNorm() const;
    Vector* toVector() const;

    void operator<<(ComplexNum& n);
    ComplexNum& operator= (const ComplexNum& n);
    virtual AbstractNum* operator+ (const AbstractNum& n) const;
    virtual AbstractNum* operator- (const AbstractNum& n) const;
    virtual AbstractNum* operator* (const AbstractNum& n) const;
    virtual AbstractNum* operator/ (const AbstractNum& n) const;
    virtual AbstractNum* operator^ (const AbstractNum& n) const;
    virtual AbstractNum* operator- () const;
    virtual bool operator== (const AbstractNum& n) const;
};

#endif
