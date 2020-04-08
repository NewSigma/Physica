/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_REALINF_H
#define PHYSICA_REALINF_H

#include "DirectedInf.h"
#include "RealNum.h"
#include "Numerical.h"
#include "Const.h"

extern const BasicConst* basicConst;

class RealInf : public DirectedInf {
public:
    bool getSign() const { return *((RealNum*)(*direction)[0])->real == basicConst->get_1(); }
    static RealInf* getInstance(bool b);

    virtual AbstractNum* operator+(const AbstractNum& n) const;
    virtual AbstractNum* operator-(const AbstractNum& n) const;
    virtual AbstractNum* operator*(const AbstractNum& n) const;
    virtual AbstractNum* operator/(const AbstractNum& n) const;
    virtual AbstractNum* operator^(const AbstractNum& n) const;
    virtual AbstractNum* operator-() const;
    virtual bool operator== (const AbstractNum& n) const;
private:
    static RealInf* positive;
    static RealInf* negative;
    RealInf(bool b);
    ~RealInf();
};

#endif
