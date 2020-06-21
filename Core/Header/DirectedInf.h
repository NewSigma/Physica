/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_DIRECTEDINF_H
#define PHYSICA_DIRECTEDINF_H

#include "AbstractNum.h"
#include "Vector.h"

namespace Physica::Core {
    class DirectedInf : public AbstractNum {
    public:
        Vector direction;

        explicit DirectedInf(const Scalar& arg);
        explicit DirectedInf(const Vector& direction);
        DirectedInf(const DirectedInf& instance);

        NumberType getType() const noexcept override;
        AbstractNum* operator+ (const AbstractNum& n) const override;
        AbstractNum* operator- (const AbstractNum& n) const override;
        AbstractNum* operator* (const AbstractNum& n) const override;
        AbstractNum* operator/ (const AbstractNum& n) const override;
        AbstractNum* operator^ (const AbstractNum& n) const override;
        AbstractNum* operator- () const override;
        bool operator== (const AbstractNum& n) const override;
    };
}



#endif
