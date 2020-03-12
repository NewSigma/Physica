/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "../../Header/DirectedInf.h"
#include "../../Header/Indeterminate.h"
#include "../../Header/ComplexInf.h"
#include "../../Header/Const.h"
#include "../../Header/RealNum.h"

extern const Const_1* const_1;

DirectedInf::DirectedInf(Vector* vec) {
    type = DirectedInfinity;
    direction = vec;
}

DirectedInf::DirectedInf(AbstractNum* arg) : DirectedInf(new Vector(arg)) {}

DirectedInf::DirectedInf(const DirectedInf& instance) : DirectedInf(instance.direction) {}

DirectedInf::DirectedInf(const DirectedInf* instance) : DirectedInf(*instance) {}

DirectedInf::~DirectedInf() {
    delete direction;
}

AbstractNum* DirectedInf::operator+ (const AbstractNum& n) const {
    switch(n.getType()) {
        case ComplexNumber:
        case RealNumber:
            return new DirectedInf(direction);
        case DirectedInfinity:
        case RealInfinity:
            if(*direction == *((DirectedInf&)n).direction)
                return new DirectedInf(this);
        case ComplexInfinity:
            return ComplexInf::getInstance();
        case Indeterminate:
            return Indeterminate::getInstance();
    }
}

AbstractNum* DirectedInf::operator- (const AbstractNum& n) const {
    switch(n.getType()) {
        case ComplexNumber:
        case RealNumber:
            return new DirectedInf(direction);
        case DirectedInfinity:
        case RealInfinity: {
            auto temp = -*direction;
            bool equal = *temp == *((DirectedInf&)n).direction;
            delete temp;
            if(equal)
                return new DirectedInf(this);
        }
        case ComplexInfinity:
            return ComplexInf::getInstance();
        case Indeterminate:
            return Indeterminate::getInstance();
    }
}

AbstractNum* DirectedInf::operator* (const AbstractNum& n) const {
    switch(n.getType()) {
        case ComplexNumber:
            return n * *this;
        case RealNumber:
        case RealInfinity:
            return new DirectedInf(direction);
        case DirectedInfinity: {
            auto arg1 = direction->toArg(0);
            auto arg2 = ((DirectedInf&)n).direction->toArg(0);
            auto result_arg = *arg1 + *arg2;
            auto result = new DirectedInf(result_arg);
            delete arg1;
            delete arg2;
            delete result_arg;
            return result;
        }
        case ComplexInfinity:
            return ComplexInf::getInstance();
        case Indeterminate:
            return Indeterminate::getInstance();
    }
}

AbstractNum* DirectedInf::operator/ (const AbstractNum& n) const {

}

AbstractNum* DirectedInf::operator^ (const AbstractNum& n) const {

}

AbstractNum* DirectedInf::operator- () const {
    return new DirectedInf(-*direction);
}

bool DirectedInf::operator== (const AbstractNum& n) const {
    if(n.getType() == DirectedInfinity) {
        direction->toUnit();
        ((DirectedInf&)n).direction->toUnit();
        return *direction == *((DirectedInf&)n).direction;
    }
    return false;
}

