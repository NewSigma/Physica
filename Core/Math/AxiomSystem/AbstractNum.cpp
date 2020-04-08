/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "../../Header/AbstractNum.h"
#include "../../Header/ComplexNum.h"
#include "../../Header/RealNum.h"
#include "../../Header/DirectedInf.h"
#include "../../Header/ComplexInf.h"
#include "../../Header/RealInf.h"
#include "../../Header/Indeterminate.h"
#include "../../Header/Numerical.h"

extern const MathConst* mathConst;

AbstractNum* AbstractNum::concretize() {
    switch(type) {
        case ComplexNumber:
            return new ComplexNum((ComplexNum*)this);
        case RealNumber:
            return new RealNum((RealNum*)this);
        case DirectedInfinity:
            return new DirectedInf((DirectedInf*)this);
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
            if(((ComplexNum&)n).imagine->length > 0)
                os << '+';
            os << *((ComplexNum&)n).imagine;
            break;
        case AbstractNum::RealNumber:
            os << *((RealNum&)n).real;
            break;
        case AbstractNum::DirectedInfinity:
            os << "DirectedInfinity" << *((DirectedInf&)n).direction;
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
    return mathConst->get_1() / n;
}

AbstractNum* sqrt(const AbstractNum& n) {
    switch(n.getType()) {
        case AbstractNum::ComplexNumber: {
            auto vec = ((ComplexNum&)n).toVector();
            auto norm = vec->toNorm();
            *((RealNum*)norm)->real << *sqrt(*((RealNum*)norm)->real);

            auto arg = vec->toArg(0);
            *((RealNum*)arg)->real /= basicConst->get_2();
            auto new_real = cos(*((RealNum*)arg)->real);
            *new_real *= *((RealNum*)norm)->real;
            auto new_imagine = sin(*((RealNum*)arg)->real);
            *new_imagine *= *((RealNum*)norm)->real;
            delete vec;
            delete norm;
            delete arg;
            return new ComplexNum(new_real, new_imagine);
        }
        case AbstractNum::RealNumber:
            return new RealNum(sqrt(*((RealNum&)n).real));
        case AbstractNum::ComplexInfinity:
            return ComplexInf::getInstance();
        case AbstractNum::RealInfinity:
            return RealInf::getInstance(((RealInf&)n).getSign());
        case AbstractNum::DirectedInfinity:
            if(((DirectedInf&)n).direction->getLength() == 2) {
                auto arg = ((DirectedInf&)n).direction->toArg(0);
                *((RealNum*)arg)->real /= basicConst->get_2();
                auto unit_vec_x = new RealNum(cos(*((RealNum*)arg)->real));
                auto unit_vec_y = new RealNum(sin(*((RealNum*)arg)->real));
                auto result = new DirectedInf(new Vector(unit_vec_x, unit_vec_y));
                delete arg;
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