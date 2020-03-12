#ifndef PHYSICA_ELEMENTARYFUNCTION_H
#define PHYSICA_ELEMENTARYFUNCTION_H

#include "AbstractNum.h"

AbstractNum* reciprocal(const AbstractNum& n);
AbstractNum* sqrt(const AbstractNum& n);
AbstractNum* factorial(const AbstractNum& n);
AbstractNum* cos(const AbstractNum& n);
AbstractNum* sin(const AbstractNum& n);
AbstractNum* tan(const AbstractNum& n);
AbstractNum* sec(const AbstractNum& n);
AbstractNum* csc(const AbstractNum& n);
AbstractNum* cot(const AbstractNum& n);
AbstractNum* arccos(const AbstractNum& n);
AbstractNum* arcsin(const AbstractNum& n);
AbstractNum* arctan(const AbstractNum& n);
AbstractNum* arcsec(const AbstractNum& n);
AbstractNum* arccsc(const AbstractNum& n);
AbstractNum* arccot(const AbstractNum& n);
AbstractNum* ln(const AbstractNum& n);
AbstractNum* log(const AbstractNum& n, const AbstractNum& a);
AbstractNum* exp(const AbstractNum& n);
AbstractNum* pow(const AbstractNum& n, const AbstractNum& a);

#endif
