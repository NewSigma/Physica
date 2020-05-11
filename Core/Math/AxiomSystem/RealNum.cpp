#include "RealNum.h"
#include "Numerical.h"
#include "DirectedInf.h"
#include "ComplexInf.h"
#include "RealInf.h"
#include "Indeterminate.h"

namespace Physica::Core {
    RealNum::RealNum(const Numerical& n) : ComplexNum(n, BasicConst::getInstance().get_0()) {}

    RealNum::RealNum(const RealNum& n) : RealNum(*n.real) {}

    RealNum::RealNum(const RealNum* n) : RealNum(*n) {}

    RealNum::~RealNum() {
        imagine = nullptr;
    }

    AbstractNum::NumberType RealNum::getType() const noexcept {
        return AbstractNum::RealNumber;
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
}