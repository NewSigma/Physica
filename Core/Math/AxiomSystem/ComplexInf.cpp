/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "../../Header/ComplexInf.h"
#include "../../Header/Indeterminate.h"

ComplexInf* ComplexInf::instance = nullptr;

ComplexInf::ComplexInf() {
    type = AbstractNum::ComplexInfinity;
    direction = nullptr;
}

ComplexInf* ComplexInf::getInstance() {
    if(instance == nullptr)
        instance = new ComplexInf();
    return instance;
}

AbstractNum* ComplexInf::operator+(const AbstractNum& n) const {
    if(n.getType() == Indeterminate)
        return Indeterminate::getInstance();
    return getInstance();
}

AbstractNum* ComplexInf::operator-(const AbstractNum& n) const {
    if(n.getType() == Indeterminate)
        return Indeterminate::getInstance();
    return getInstance();
}

AbstractNum* ComplexInf::operator*(const AbstractNum& n) const {
    if(n.getType() == Indeterminate)
        return Indeterminate::getInstance();
    return getInstance();
}

AbstractNum* ComplexInf::operator/(const AbstractNum& n) const {
    if(n.getType() == ComplexNumber || n.getType() == RealNumber)
        return getInstance();
    return Indeterminate::getInstance();
}

AbstractNum* ComplexInf::operator^(const AbstractNum& n) const {
    if(n.getType() == ComplexNumber || n.getType() == RealNumber)
        return getInstance();
    return Indeterminate::getInstance();
}

AbstractNum* ComplexInf::operator-() const {
    return instance;
}

bool ComplexInf::operator== (const AbstractNum& n) const {
    return n.getType() == ComplexInfinity;
}