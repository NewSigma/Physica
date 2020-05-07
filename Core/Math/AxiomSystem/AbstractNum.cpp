/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "AbstractNum.h"
#include "ComplexNum.h"
#include "RealNum.h"
#include "DirectedInf.h"
#include "ComplexInf.h"
#include "RealInf.h"
#include "Indeterminate.h"
#include "Numerical.h"
#include "ElementaryFunction.h"

AbstractNum* AbstractNum::concretize() {
    switch(getType()) {
        case ComplexNumber:
            return new ComplexNum((ComplexNum*)this);
        case RealNumber:
            return new RealNum((RealNum*)this);
        case DirectedInfinity:
            return new DirectedInf(*(DirectedInf*)this);
        case ComplexInfinity:
            return ComplexInf::getInstance();
        case RealInfinity:
            return RealInf::getInstance(((RealInf*)this)->getSign());
        case Indeterminate:
            return Indeterminate::getInstance();
    }
}

std::ostream& operator<<(std::ostream& os, const AbstractNum& n) {
    switch(n.getType()) {
        case AbstractNum::ComplexNumber:
            os << ((ComplexNum&)n).real;
            if(((ComplexNum&)n).imagine->getLength() > 0)
                os << '+';
            os << *((ComplexNum&)n).imagine;
            break;
        case AbstractNum::RealNumber:
            os << *((RealNum&)n).real;
            break;
        case AbstractNum::DirectedInfinity:
            os << "DirectedInfinity" << ((DirectedInf&)n).direction;
            break;
        case AbstractNum::ComplexInfinity:
            os << "ComplexInfinity";
            break;
        case AbstractNum::RealInfinity:
            if(((RealInf&)n).getSign())
                os << "Infinity";
            else
                os << "-Infinity";
            break;
        case AbstractNum::Indeterminate:
            os << "Indeterminate";
            break;
    }
    return os;
}

bool AbstractNum::operator!= (const AbstractNum& n) const {
    return !(*this == n);
}
////////////////////////////////////////Elementary Functions////////////////////////////////////////////
AbstractNum* reciprocal(const AbstractNum& n) {
    return MathConst::getInstance().get_1() / n;
}

AbstractNum* sqrt(const AbstractNum& n) {
    switch(n.getType()) {
        case AbstractNum::ComplexNumber: {
            Vector vec = ((ComplexNum&)n).toVector();
            Numerical norm = sqrt(vec.toNorm());
            Numerical arg = vec.toArg(0);
            return new ComplexNum(norm * cos(arg), norm * sin(arg));
        }
        case AbstractNum::RealNumber:
            return new RealNum(sqrt(*((RealNum&)n).real));
        case AbstractNum::ComplexInfinity:
            return ComplexInf::getInstance();
        case AbstractNum::RealInfinity:
            return RealInf::getInstance(((RealInf&)n).getSign());
        case AbstractNum::DirectedInfinity:
            if(((DirectedInf&)n).direction.getLength() == 2) {
                Numerical arg = ((DirectedInf&)n).direction.toArg(0);
                arg /= BasicConst::getInstance().get_2();
                auto result = new DirectedInf(Vector(cos(arg), sin(arg)));
                return result;
            }
        case AbstractNum::Indeterminate:
            return Indeterminate::getInstance();
    }
}

AbstractNum* factorial(const AbstractNum& n) {

}

AbstractNum* ln(const AbstractNum& n) {

}

AbstractNum* log(const AbstractNum& n, const AbstractNum& a) {

}

AbstractNum* exp(const AbstractNum& n) {

}

AbstractNum* pow(const AbstractNum& n, const AbstractNum& a) {

}
AbstractNum* cos(const AbstractNum& n) {

}

AbstractNum* sin(const AbstractNum& n) {

}

AbstractNum* tan(const AbstractNum& n) {

}

AbstractNum* sec(const AbstractNum& n) {

}

AbstractNum* csc(const AbstractNum& n) {

}

AbstractNum* cot(const AbstractNum& n) {

}

AbstractNum* arccos(const AbstractNum& n) {

}

AbstractNum* arcsin(const AbstractNum& n) {

}

AbstractNum* arctan(const AbstractNum& n) {

}

AbstractNum* arcsec(const AbstractNum& n) {

}

AbstractNum* arccsc(const AbstractNum& n) {

}

AbstractNum* arccot(const AbstractNum& n) {

}

AbstractNum* cosh(const AbstractNum& n) {

}

AbstractNum* sinh(const AbstractNum& n) {

}

AbstractNum* tanh(const AbstractNum& n) {

}

AbstractNum* sech(const AbstractNum& n) {

}

AbstractNum* csch(const AbstractNum& n) {

}

AbstractNum* coth(const AbstractNum& n) {

}

AbstractNum* arccosh(const AbstractNum& n) {

}

AbstractNum* arcsinh(const AbstractNum& n) {

}

AbstractNum* arctanh(const AbstractNum& n) {

}

AbstractNum* arcsech(const AbstractNum& n) {

}

AbstractNum* arccsch(const AbstractNum& n) {

}

AbstractNum* arccoth(const AbstractNum& n) {

}