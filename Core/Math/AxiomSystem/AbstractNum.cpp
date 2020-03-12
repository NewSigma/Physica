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
            if(((ComplexNum&)n).imagine->sign)
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