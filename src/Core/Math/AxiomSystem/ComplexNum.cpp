/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "Physica/Core/ComplexNum.h"
#include "Physica/Core/RealNum.h"
#include "Physica/Core/DirectedInf.h"
#include "Physica/Core/ComplexInf.h"
#include "Physica/Core/RealInf.h"
#include "Physica/Core/Indeterminate.h"
#include "Physica/Core/ElementaryFunction.h"

namespace Physica::Core {
    ComplexNum::ComplexNum(const Scalar& n1, const Scalar& n2, bool polar) {
        if(polar) {
            real = new Scalar(n1 * cos(n2));
            imagine = new Scalar(n1 * sin(n2));
        }
        else {
            real = new Scalar(n1);
            imagine = new Scalar(n2);
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

    Scalar ComplexNum::toNorm() const {
        return *real * *real + *imagine * *imagine;
    }

    Vector ComplexNum::toVector() const {
        Vector result(2);
        result.grow(*real);
        result.grow(*imagine);
        return result;
    }

    AbstractNum::NumberType ComplexNum::getType() const noexcept {
        return AbstractNum::ComplexNumber;
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
        real = new Scalar(*n.real);
        imagine = new Scalar(*n.imagine);
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
                Vector unit = toVector();
                auto arg1 = unit.toArg(0);
                auto arg2 = ((DirectedInf&)n).direction.toArg(0);
                auto result_arg = arg1 + arg2;
                return new DirectedInf(result_arg);
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
                Scalar new_real = *real * *((ComplexNum&)n).real + *imagine * *((ComplexNum&)n).imagine;
                Scalar new_imagine = *((ComplexNum&)n).real * *imagine - *real * *((ComplexNum&)n).imagine;
                Scalar temp = *((ComplexNum&)n).imagine * *((ComplexNum&)n).imagine + *((ComplexNum&)n).real * *((ComplexNum&)n).real;
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
}