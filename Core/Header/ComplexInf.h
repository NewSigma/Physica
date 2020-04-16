/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_COMPLEXINF_H
#define PHYSICA_COMPLEXINF_H

#include "DirectedInf.h"

class ComplexInf : public DirectedInf {
public:
    static ComplexInf* getInstance();

    NumberType getType() const noexcept override;
    AbstractNum* operator+(const AbstractNum& n) const override;
    AbstractNum* operator-(const AbstractNum& n) const override;
    AbstractNum* operator*(const AbstractNum& n) const override;
    AbstractNum* operator/(const AbstractNum& n) const override;
    AbstractNum* operator^(const AbstractNum& n) const override;
    AbstractNum* operator-() const override;
    bool operator== (const AbstractNum& n) const override;
private:
    static ComplexInf* instance;
    ComplexInf();
};

#endif
