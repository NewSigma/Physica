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

    ComplexNum(const Numerical& n1, const Numerical& n2, bool polar = false);
    ComplexNum(ComplexNum& instance);
    explicit ComplexNum(ComplexNum* instance);
    ~ComplexNum() override;
    ComplexNum* toConjugate() const;
    Numerical toNorm() const;
    Vector toVector() const;

    NumberType getType() const noexcept override;
    void operator<<(ComplexNum& n);
    ComplexNum& operator= (const ComplexNum& n);
    AbstractNum* operator+ (const AbstractNum& n) const override;
    AbstractNum* operator- (const AbstractNum& n) const override;
    AbstractNum* operator* (const AbstractNum& n) const override;
    AbstractNum* operator/ (const AbstractNum& n) const override;
    AbstractNum* operator^ (const AbstractNum& n) const override;
    AbstractNum* operator- () const override;
    bool operator== (const AbstractNum& n) const override;
};

#endif
