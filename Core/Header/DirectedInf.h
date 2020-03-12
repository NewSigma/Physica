/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_DIRECTEDINF_H
#define PHYSICA_DIRECTEDINF_H

#include "AbstractNum.h"
#include "Vector.h"

class DirectedInf : public AbstractNum {
public:
    Vector* direction;

    DirectedInf(Vector* direction);
    DirectedInf(AbstractNum* arg);
    DirectedInf(const DirectedInf& instance);
    DirectedInf(const DirectedInf* instance);
    ~DirectedInf();

    virtual AbstractNum* operator+ (const AbstractNum& n) const;
    virtual AbstractNum* operator- (const AbstractNum& n) const;
    virtual AbstractNum* operator* (const AbstractNum& n) const;
    virtual AbstractNum* operator/ (const AbstractNum& n) const;
    virtual AbstractNum* operator^ (const AbstractNum& n) const;
    virtual AbstractNum* operator- () const;
    virtual bool operator== (const AbstractNum& n) const;
protected:
    DirectedInf() = default;
};

#endif
