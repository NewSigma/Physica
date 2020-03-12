/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_INDETERMINATE_H
#define PHYSICA_INDETERMINATE_H

#include "AbstractNum.h"

class Indeterminate : public AbstractNum {
public:
    static Indeterminate* getInstance();
    virtual AbstractNum* operator+(const AbstractNum& n) const;
    virtual AbstractNum* operator-(const AbstractNum& n) const;
    virtual AbstractNum* operator*(const AbstractNum& n) const;
    virtual AbstractNum* operator/(const AbstractNum& n) const;
    virtual AbstractNum* operator^(const AbstractNum& n) const;
    virtual AbstractNum* operator-() const;
    virtual bool operator== (const AbstractNum& n) const;
private:
    static Indeterminate* instance;
    Indeterminate();
    virtual ~Indeterminate();
};

#endif
