/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <Physica/Core/Scalar.h>
#include "Physica/Core/DirectedInf.h"
#include "Physica/Core/Indeterminate.h"
#include "Physica/Core/ComplexInf.h"
#include "Physica/Core/Vector.h"

namespace Physica::Core {
    DirectedInf::DirectedInf(const Scalar& arg) : direction(2) {
        direction.grow(cos(arg));
        direction.grow(sin(arg));
    }

    DirectedInf::DirectedInf(const Vector& vec) : direction(Vector(vec)) {
        direction.toUnit();
    }

    DirectedInf::DirectedInf(const DirectedInf& instance) : DirectedInf(Vector(instance.direction)) {}

    AbstractNum::NumberType DirectedInf::getType() const noexcept {
        return AbstractNum::DirectedInfinity;
    }

    AbstractNum* DirectedInf::operator+ (const AbstractNum& n) const {
        switch(n.getType()) {
            case ComplexNumber:
            case RealNumber:
                return new DirectedInf(direction);
            case DirectedInfinity:
            case RealInfinity:
                if(direction == ((DirectedInf&)n).direction)
                    return new DirectedInf(*this);
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
                Vector temp = -direction;
                bool equal = temp == ((DirectedInf&)n).direction;
                if(equal)
                    return new DirectedInf(*this);
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
                Scalar arg1 = direction.toArg(0);
                Scalar arg2 = ((DirectedInf&)n).direction.toArg(0);
                Scalar result_arg = arg1 + arg2;
                auto result = new DirectedInf(result_arg);
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
        return new DirectedInf(-direction);
    }

    bool DirectedInf::operator== (const AbstractNum& n) const {
        if(n.getType() == DirectedInfinity)
            return direction == ((DirectedInf&)n).direction;
        return false;
    }
}