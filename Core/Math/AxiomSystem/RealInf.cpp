/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "../../Header/RealInf.h"
#include "../../Header/ComplexInf.h"
#include "../../Header/Indeterminate.h"

RealInf* RealInf::positive = nullptr;
RealInf* RealInf::negative = nullptr;

RealInf::RealInf(bool b) {
    type = AbstractNum::RealInfinity;
    auto arr = new AbstractNum*[1];
    if(b)
        arr[0] = new RealNum(new Numerical(basicConst->get_1()));
    else
        arr[0] = new RealNum(new Numerical(basicConst->getMinus_1()));
    direction = new Vector(arr, 1);
}

RealInf::~RealInf() {
    direction->numbers[0] = nullptr;
    delete direction;
}

RealInf* RealInf::getInstance(bool b) {
    if(b) {
        if(positive == nullptr)
            positive = new RealInf(b);
        return positive;
    }
    else {
        if(negative == nullptr)
            negative = new RealInf(b);
        return negative;
    }
}

AbstractNum* RealInf::operator+(const AbstractNum& n) const {
    switch(n.getType()) {
        case AbstractNum::ComplexNumber:
        case AbstractNum::RealNumber:
            return (AbstractNum*)this;
        case AbstractNum::DirectedInfinity:
            if((*((DirectedInf&)n).direction) == *direction)
                return (AbstractNum*)this;
            return ComplexInf::getInstance();
        case AbstractNum::RealInfinity:
            if(((RealInf&)n).getSign() == getSign())
                return (AbstractNum*)this;
        case AbstractNum::ComplexInfinity:
        case AbstractNum::Indeterminate:
            return Indeterminate::getInstance();
    }
}

AbstractNum* RealInf::operator-(const AbstractNum& n) const {
    switch(n.getType()) {
        case AbstractNum::ComplexNumber:
        case AbstractNum::RealNumber:
            return (AbstractNum*)this;
        case AbstractNum::DirectedInfinity:
            if((*((DirectedInf&)n).direction) == *direction)
                return (AbstractNum*)this;
            return ComplexInf::getInstance();
        case AbstractNum::RealInfinity:
            if(((RealInf&)n).getSign() != getSign())
                return (AbstractNum*)this;
        case AbstractNum::ComplexInfinity:
        case AbstractNum::Indeterminate:
            return Indeterminate::getInstance();
    }
}

AbstractNum* RealInf::operator*(const AbstractNum& n) const {
    switch(n.getType()) {
        case ComplexNumber:
            return new DirectedInf(((ComplexNum&)n).toVector());
        case DirectedInfinity:
            return new DirectedInf((DirectedInf&)n);
        case RealNumber:
            return getInstance(getSign() && ((RealNum&)n).real->getLength() > 0);
        case RealInfinity:
            return getInstance(getSign() && ((RealInf&)n).getSign());
        case ComplexInfinity:
            return ComplexInf::getInstance();
        case Indeterminate:
            return Indeterminate::getInstance();
    }
}

AbstractNum* RealInf::operator/(const AbstractNum& n) const {
    switch(n.getType()) {
        case ComplexNumber: {
            auto temp = reciprocal(n);
            auto result = new DirectedInf(((ComplexNum*)temp)->toVector());
            delete temp;
            return result;
        }
        case RealNumber:
            return getInstance(getSign() && ((RealNum&)n).real->getLength() > 0);
        case DirectedInfinity:
        case ComplexInfinity:
        case RealInfinity:
        case Indeterminate:
            return Indeterminate::getInstance();
    }
}

AbstractNum* RealInf::operator^(const AbstractNum& n) const {
    switch(n.getType()) {
        case RealNumber:
            if(((RealNum&)n).real->getLength() > 0)
                return getInstance(getSign() && ((RealNum&)n).real->getLength() > 0);
            else
                return new RealNum(getZero());
        case ComplexNumber:
        case DirectedInfinity:
        case RealInfinity:
            return ComplexInf::getInstance();
        case ComplexInfinity:
        case Indeterminate:
            return Indeterminate::getInstance();
    }
}

AbstractNum* RealInf::operator-() const {
    return getInstance(!getSign());
}

bool RealInf::operator== (const AbstractNum& n) const {
    if(n.getType() == RealInfinity) {
        return getSign() == ((RealInf&)n).getSign();
    }
    return false;
}