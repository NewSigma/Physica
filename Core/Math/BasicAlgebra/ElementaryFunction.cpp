/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "../../Header/Numerical.h"
#include "../../Header/ElementaryFunction.h"
#include "../../Header/Indeterminate.h"
#include "../../Header/ComplexInf.h"
#include "../../Header/RealInf.h"

extern const Const_2* const_2;

AbstractNum* reciprocal(const AbstractNum& n) {
    return *const_2->_1 / n;
}

AbstractNum* sqrt(const AbstractNum& n) {
    switch(n.getType()) {
        case AbstractNum::ComplexNumber: {
            auto vec = ((ComplexNum&)n).toVector();
            auto norm = vec->toNorm();
            *((RealNum*)norm)->real << *sqrt(*((RealNum*)norm)->real);

            auto arg = vec->toArg(0);
            *((RealNum*)arg)->real /= *const_1->_2;
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
                *((RealNum*)arg)->real /= *const_1->_2;
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

AbstractNum* ln(const AbstractNum& n) {

}

AbstractNum* log(const AbstractNum& n, const AbstractNum& a) {

}

AbstractNum* exp(const AbstractNum& n) {

}

AbstractNum* pow(const AbstractNum& n, const AbstractNum& a) {

}