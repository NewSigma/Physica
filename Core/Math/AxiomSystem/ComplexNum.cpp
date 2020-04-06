/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "../../Header/ComplexNum.h"
#include "../../Header/RealNum.h"
#include "../../Header/DirectedInf.h"
#include "../../Header/ComplexInf.h"
#include "../../Header/RealInf.h"
#include "../../Header/Indeterminate.h"
#include "../../Header/Vector.h"
#include "../../Header/Numerical.h"

ComplexNum::ComplexNum(Numerical* n1, Numerical* n2, bool polar) {
    type = AbstractNum::ComplexNumber;
    if(polar) {
        real = cos(*n2);
        *real *= *n1;
        imagine = sin(*n2);
        *imagine *= *n1;
        delete n2;
    }
    else {
        real = n1;
        imagine = n2;
    }
}

ComplexNum::ComplexNum(ComplexNum& instance) : ComplexNum(new Numerical(instance.real), new Numerical(instance.imagine)) {}

ComplexNum::ComplexNum(ComplexNum* instance) : ComplexNum(*instance) {}

ComplexNum::~ComplexNum() {
    delete real;
    delete imagine;
}

ComplexNum* ComplexNum::toConjugate() const {
    return new ComplexNum(new Numerical(real), -*imagine);
}

Numerical* ComplexNum::toNorm() const {
    auto result = *real * *real;
    auto temp = *imagine * *imagine;
    *result += *temp;
    delete temp;
    return result;
}

Vector* ComplexNum::toVector() const {
    auto arr = new AbstractNum*[2];
    arr[0] = new RealNum(new Numerical(real));
    arr[1] = new RealNum(new Numerical(imagine));
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
            return new ComplexNum(*real + *((RealNum&)n).real, new Numerical(imagine));
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
            auto new_real = *real * *((ComplexNum&)n).real;
            auto new_imagine = *real * *((ComplexNum&)n).imagine;

            auto temp = *imagine * *((ComplexNum&)n).imagine;
            *new_real -= *temp;
            delete temp;

            temp = *imagine * *((ComplexNum&)n).real;
            *new_imagine += *temp;
            delete temp;
            return new ComplexNum(new_real, new_imagine);
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
            auto new_real = *real * *((ComplexNum&)n).real;
            auto new_imagine = *((ComplexNum&)n).real * *imagine;

            auto temp = *imagine * *((ComplexNum&)n).imagine;
            *new_real += *temp;
            delete temp;

            temp = *real * *((ComplexNum&)n).imagine;
            *new_imagine -= *temp;
            delete temp;

            auto temp1 = *((ComplexNum&)n).imagine * *((ComplexNum&)n).imagine;
            temp = *((ComplexNum&)n).real * *((ComplexNum&)n).real;
            *temp += *temp1;
            delete temp1;

            *new_real /= *temp;
            *new_imagine /= *temp;
            delete temp;
            return new ComplexNum(new_real, new_imagine);
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