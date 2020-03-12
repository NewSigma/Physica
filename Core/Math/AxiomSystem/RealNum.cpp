#include "../../Header/RealNum.h"
#include "../../Header/Numerical.h"
#include "../../Header/DirectedInf.h"
#include "../../Header/ComplexInf.h"
#include "../../Header/RealInf.h"
#include "../../Header/Indeterminate.h"

RealNum::RealNum(Numerical* n) : ComplexNum(n, (Numerical*)const_1->_0) {
    type = RealNumber;
}

RealNum::RealNum(const RealNum& n) : RealNum(new Numerical(n.real)) {}

RealNum::RealNum(const RealNum* n) : RealNum(*n) {}

RealNum::~RealNum() {
    imagine = nullptr;
}

AbstractNum* RealNum::operator+ (const AbstractNum& n) const {
    switch(n.getType()) {
        case ComplexNumber:
            return ((ComplexNum&)n) + *this;
        case RealNumber:
            return new RealNum(*real + *((RealNum&)n).real);
        case DirectedInfinity:
            return ((DirectedInf&)n) + *this;
        case ComplexInfinity:
            return ((ComplexInf&)n) + *this;
        case RealInfinity:
            return ((RealInf&)n) + *this;
        case Indeterminate:
            return Indeterminate::getInstance();
    }
}

AbstractNum* RealNum::operator- (const AbstractNum& n) const {
    switch(n.getType()) {
        case ComplexNumber:
            return ((ComplexNum&)n) - *this;
        case RealNumber:
            return new RealNum(*real - *((RealNum&)n).real);
        case DirectedInfinity:
            return ((DirectedInf&)n) - *this;
        case ComplexInfinity:
            return ((ComplexInf&)n) - *this;
        case RealInfinity:
            return ((RealInf&)n) - *this;
        case Indeterminate:
            return Indeterminate::getInstance();
    }
}

AbstractNum* RealNum::operator* (const AbstractNum& n) const {
    switch(n.getType()) {
        case ComplexNumber:
            return ((ComplexNum&)n) * *this;
        case RealNumber:
            return new RealNum(*real * *((RealNum&)n).real);
        case DirectedInfinity:
            return ((DirectedInf&)n) * *this;
        case ComplexInfinity:
            return ((ComplexInf&)n) * *this;
        case RealInfinity:
            return ((RealInf&)n) * *this;
        case Indeterminate:
            return Indeterminate::getInstance();
    }
}

AbstractNum* RealNum::operator/ (const AbstractNum& n) const {
    switch(n.getType()) {
        case ComplexNumber:
            return ((ComplexNum&)n) / *this;
        case RealNumber:
            return new RealNum(*real / *((RealNum&)n).real);
        case DirectedInfinity:
            return ((DirectedInf&)n) / *this;
        case ComplexInfinity:
            return ((ComplexInf&)n) / *this;
        case RealInfinity:
            return ((RealInf&)n) / *this;
        case Indeterminate:
            return Indeterminate::getInstance();
    }
}

AbstractNum* RealNum::operator^ (const AbstractNum& n) const {
    switch(n.getType()) {
        case ComplexNumber:
            return ((ComplexNum&)n) ^ *this;
        case RealNumber:
            return new RealNum(*real ^ *((RealNum&)n).real);
        case DirectedInfinity:
            return ((DirectedInf&)n) ^ *this;
        case ComplexInfinity:
            return ((ComplexInf&)n) ^ *this;
        case RealInfinity:
            return ((RealInf&)n) ^ *this;
        case Indeterminate:
            return Indeterminate::getInstance();
    }
}

AbstractNum* RealNum::operator- () const {
    return new RealNum(-*real);
}

bool RealNum::operator== (const AbstractNum& n) const {
    if(n.getType() == RealNumber) {
        return *real == *((RealNum&)n).real;
    }
    return false;
}

bool RealNum::operator> (const RealNum& n) const {
    return *real > *n.real;
}

bool RealNum::operator< (const RealNum& n) const {
    return *real < *n.real;
}

bool RealNum::operator>= (const RealNum& n) const {
    return *real >= *n.real;
}

bool RealNum::operator<= (const RealNum& n) const {
    return *real <= *n.real;
}