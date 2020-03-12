/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_ABSTRACTNUM_H
#define PHYSICA_ABSTRACTNUM_H

#include <ostream>

//AbstractNum is a abstract class.
class AbstractNum {
public:
    enum NumberType {
        ComplexNumber,
        RealNumber,
        DirectedInfinity,
        ComplexInfinity,
        RealInfinity,
        Indeterminate
    };
    virtual ~AbstractNum() = default;
    NumberType getType() const { return type; }
    AbstractNum* concretize();

    friend std::ostream& operator<<(std::ostream& os, const AbstractNum& n);
    virtual AbstractNum* operator+ (const AbstractNum& n) const = 0;
    virtual AbstractNum* operator- (const AbstractNum& n) const = 0;
    virtual AbstractNum* operator* (const AbstractNum& n) const = 0;
    virtual AbstractNum* operator/ (const AbstractNum& n) const = 0;
    virtual AbstractNum* operator^ (const AbstractNum& n) const = 0;
    virtual AbstractNum* operator- () const = 0;
    virtual bool operator== (const AbstractNum& n) const = 0;
    bool operator!= (const AbstractNum& n) const;
protected:
    NumberType type;
};

#endif
