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

    NumberType getType() const noexcept override;
    AbstractNum* operator+(const AbstractNum& n) const override;
    AbstractNum* operator-(const AbstractNum& n) const override;
    AbstractNum* operator*(const AbstractNum& n) const override;
    AbstractNum* operator/(const AbstractNum& n) const override;
    AbstractNum* operator^(const AbstractNum& n) const override;
    AbstractNum* operator-() const override;
    bool operator== (const AbstractNum& n) const override;
private:
    static RealInf* positive;
    static RealInf* negative;
    RealInf(bool b);
    ~RealInf();
};

#endif
