/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_INDETERMINATE_H
#define PHYSICA_INDETERMINATE_H

#include "AbstractNum.h"

namespace Physica::Core {
    class Indeterminate : public AbstractNum {
    public:
        static Indeterminate* getInstance();
        NumberType getType() const noexcept override;
        AbstractNum* operator+(const AbstractNum& n) const override;
        AbstractNum* operator-(const AbstractNum& n) const override;
        AbstractNum* operator*(const AbstractNum& n) const override;
        AbstractNum* operator/(const AbstractNum& n) const override;
        AbstractNum* operator^(const AbstractNum& n) const override;
        AbstractNum* operator-() const override;
        bool operator== (const AbstractNum& n) const override;
    private:
        static Indeterminate* instance;
        Indeterminate();
        ~Indeterminate();
    };
}

#endif
