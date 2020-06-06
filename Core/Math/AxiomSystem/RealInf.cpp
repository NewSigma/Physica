/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "RealInf.h"
#include "ComplexInf.h"
#include "Indeterminate.h"

namespace Physica::Core {
    RealInf* RealInf::positive = nullptr;
    RealInf* RealInf::negative = nullptr;

    RealInf::RealInf(bool b) : DirectedInf(Vector()) {
        direction = Vector(1);
        Numerical temp = b ? BasicConst::getInstance().get_1() : BasicConst::getInstance().getMinus_1();
        direction.grow(temp);
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

    AbstractNum::NumberType RealInf::getType() const noexcept {
        return AbstractNum::RealInfinity;
    }

    AbstractNum* RealInf::operator+(const AbstractNum& n) const {
        switch(n.getType()) {
            case AbstractNum::ComplexNumber:
            case AbstractNum::RealNumber:
                return (AbstractNum*)this;
            case AbstractNum::DirectedInfinity:
                if((((DirectedInf&)n).direction) == direction)
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
                if((((DirectedInf&)n).direction) == direction)
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
}