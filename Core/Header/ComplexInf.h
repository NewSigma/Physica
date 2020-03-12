/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_COMPLEXINF_H
#define PHYSICA_COMPLEXINF_H

#include "DirectedInf.h"

class ComplexInf : public DirectedInf {
public:
    static ComplexInf* getInstance();

    virtual AbstractNum* operator+(const AbstractNum& n) const;
    virtual AbstractNum* operator-(const AbstractNum& n) const;
    virtual AbstractNum* operator*(const AbstractNum& n) const;
    virtual AbstractNum* operator/(const AbstractNum& n) const;
    virtual AbstractNum* operator^(const AbstractNum& n) const;
    virtual AbstractNum* operator-() const;
    virtual bool operator== (const AbstractNum& n) const;
private:
    static ComplexInf* instance;
    ComplexInf();
};

#endif
