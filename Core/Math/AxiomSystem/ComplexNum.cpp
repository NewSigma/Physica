/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "ComplexNum.h"
#include "RealNum.h"
#include "DirectedInf.h"
#include "ComplexInf.h"
#include "RealInf.h"
#include "Indeterminate.h"
#include "ElementaryFunction.h"

ComplexNum::ComplexNum(const Numerical& n1, const Numerical& n2, bool polar) {
    type = AbstractNum::ComplexNumber;
    if(polar) {
        real = new Numerical(n1 * cos(n2));
        imagine = new Numerical(n1 * sin(n2));
    }
    else {
        real = new Numerical(n1);
        imagine = new Numerical(n2);
    }
}

ComplexNum::ComplexNum(ComplexNum& instance) : ComplexNum(*instance.real, *instance.imagine) {}

ComplexNum::ComplexNum(ComplexNum* instance) : ComplexNum(*instance) {}

ComplexNum::~ComplexNum() {
    delete real;
    delete imagine;
}

ComplexNum* ComplexNum::toConjugate() const {
    return new ComplexNum(*real, -*imagine);
}

Numerical ComplexNum::toNorm() const {
    return *real * *real + *imagine * *imagine;
}

Vector* ComplexNum::toVector() const {
    auto arr = new AbstractNum*[2];
    arr[0] = new RealNum(*real);
    arr[1] = new RealNum(*imagine);
    return new Vector(arr, 2);
}

void ComplexNum::operator<<(ComplexNum& n) {
    this->~ComplexNum();
    real = n.real;
    imagine = n.imagine;
}

ComplexNum& ComplexNum::operator= (const ComplexNum& n) {
    if(this == &n)
        return *this;
    this->~ComplexNum();
    real = new Numerical(n.real);
    imagine = new Numerical(n.imagine);
    return *this;
}

AbstractNum* ComplexNum::operator+ (const AbstractNum& n) const {
    switch(n.getType()) {
        case ComplexNumber:
            return new ComplexNum(*real + *((ComplexNum&)n).real, *imagine + *((ComplexNum&)n).imagine);
        case RealNumber:
            return new ComplexNum(*real + *((RealNum&)n).real, *imagine);
        case DirectedInfinity:
            return new DirectedInf((DirectedInf&)n);
        case ComplexInfinity:
            return ComplexInf::getInstance();
        case RealInfinity:
            return RealInf::getInstance(((RealInf&)n).getSign());
        case Indeterminate:
            return Indeterminate::getInstance();
    }
}

AbstractNum* ComplexNum::operator- (const AbstractNum& n) const {
    switch(n.getType()) {
        case ComplexNumber:
            return new ComplexNum(*real - *((ComplexNum&)n).real, *imagine - *((ComplexNum&)n).imagine);
        case RealNumber:
            return new ComplexNum(*real - *((RealNum&)n).real, -*imagine);
        case DirectedInfinity:
            return -(DirectedInf&)n;
        case ComplexInfinity:
            return ComplexInf::getInstance();
        case RealInfinity:
            return RealInf::getInstance(!((RealInf&)n).getSign());
        case Indeterminate:
            return Indeterminate::getInstance();
    }
}

AbstractNum* ComplexNum::operator* (const AbstractNum& n) const {
    switch(n.getType()) {
        case ComplexNumber: {
            return new ComplexNum(*real * *((ComplexNum&)n).real - *imagine * *((ComplexNum&)n).imagine
                    , *real * *((ComplexNum&)n).imagine + *imagine * *((ComplexNum&)n).real);
        }
        case RealNumber:
            return new ComplexNum(*real * *((RealNum&)n).real, *imagine * *((RealNum&)n).real);
        case DirectedInfinity: {
            auto unit = toVector();
            auto arg1 = unit->toArg(0);
            auto arg2 = ((DirectedInf&)n).direction->toArg(0);
            auto result_arg = *arg1 + *arg2;
            auto result = new DirectedInf(result_arg);
            delete unit;
            delete arg1;
            delete arg2;
            delete result_arg;
            return result;
        }
        case ComplexInfinity:
            return ComplexInf::getInstance();
        case RealInfinity:
            return new DirectedInf(toVector());
        case Indeterminate:
            return Indeterminate::getInstance();
    }
}

AbstractNum* ComplexNum::operator/ (const AbstractNum& n) const {
    switch(n.getType()) {
        case ComplexNumber: {
            Numerical new_real = *real * *((ComplexNum&)n).real + *imagine * *((ComplexNum&)n).imagine;
            Numerical new_imagine = *((ComplexNum&)n).real * *imagine - *real * *((ComplexNum&)n).imagine;
            Numerical temp = *((ComplexNum&)n).imagine * *((ComplexNum&)n).imagine + *((ComplexNum&)n).real * *((ComplexNum&)n).real;
            return new ComplexNum(new_real / temp, new_imagine / temp);
        }
        case RealNumber:
            return new ComplexNum(*real / *((RealNum&)n).real, *imagine / *((RealNum&)n).real);
        case DirectedInfinity:
        case ComplexInfinity:
        case RealInfinity:
            return new RealNum(getZero());
        case Indeterminate:
            return Indeterminate::getInstance();
    }
}

AbstractNum* ComplexNum::operator^ (const AbstractNum& n) const {

}

AbstractNum* ComplexNum::operator- () const {
    return new ComplexNum(-*real, -*imagine);
}

bool ComplexNum::operator== (const AbstractNum& n) const {
    if(n.getType() == ComplexNumber)
        return *real == *((ComplexNum&)n).real && *imagine == *((ComplexNum&)n).imagine;
    return false;
}