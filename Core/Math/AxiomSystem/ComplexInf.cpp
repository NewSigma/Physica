/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "Core/Header/ComplexInf.h"
#include "Core/Header/Indeterminate.h"

namespace Physica::Core {
    ComplexInf* ComplexInf::instance = nullptr;

    ComplexInf::ComplexInf() : DirectedInf(Vector()) {}

    ComplexInf* ComplexInf::getInstance() {
        if(instance == nullptr)
            instance = new ComplexInf();
        return instance;
    }

    AbstractNum::NumberType ComplexInf::getType() const noexcept {
        return AbstractNum::ComplexInfinity;
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
}