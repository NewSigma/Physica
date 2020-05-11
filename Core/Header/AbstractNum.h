/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_ABSTRACTNUM_H
#define PHYSICA_ABSTRACTNUM_H

#include <ostream>

namespace Physica::Core {
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
        virtual NumberType getType() const noexcept = 0;
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
    };
////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    AbstractNum* reciprocal(const AbstractNum& n);
    AbstractNum* sqrt(const AbstractNum& n);
    AbstractNum* factorial(const AbstractNum& n);
    AbstractNum* ln(const AbstractNum& n);
    AbstractNum* log(const AbstractNum& n, const AbstractNum& a);
    AbstractNum* exp(const AbstractNum& n);
    AbstractNum* pow(const AbstractNum& n, const AbstractNum& a);
    AbstractNum* cos(const AbstractNum& n);
    AbstractNum* sin(const AbstractNum& n);
    AbstractNum* tan(const AbstractNum& n);
    AbstractNum* sec(const AbstractNum& n);
    AbstractNum* csc(const AbstractNum& n);
    AbstractNum* cot(const AbstractNum& n);
    AbstractNum* arccos(const AbstractNum& n);
    AbstractNum* arcsin(const AbstractNum& n);
    AbstractNum* arctan(const AbstractNum& n);
    AbstractNum* arcsec(const AbstractNum& n);
    AbstractNum* arccsc(const AbstractNum& n);
    AbstractNum* arccot(const AbstractNum& n);
    AbstractNum* cosh(const AbstractNum& n);
    AbstractNum* sinh(const AbstractNum& n);
    AbstractNum* tanh(const AbstractNum& n);
    AbstractNum* sech(const AbstractNum& n);
    AbstractNum* csch(const AbstractNum& n);
    AbstractNum* coth(const AbstractNum& n);
    AbstractNum* arccosh(const AbstractNum& n);
    AbstractNum* arcsinh(const AbstractNum& n);
    AbstractNum* arctanh(const AbstractNum& n);
    AbstractNum* arcsech(const AbstractNum& n);
    AbstractNum* arccsch(const AbstractNum& n);
    AbstractNum* arccoth(const AbstractNum& n);
}

#endif
