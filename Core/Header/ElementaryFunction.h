/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_ELEMENTARYFUNCTION_H
#define PHYSICA_ELEMENTARYFUNCTION_H

class Numerical;

Numerical* randomNumerical();
Numerical* randomNumerical(Numerical* lowerBound, Numerical* upperBound);
Numerical* reciprocal(const Numerical& n);
Numerical* sqrt_light(const Numerical& n);
Numerical* sqrt(const Numerical& n);
Numerical* factorial(const Numerical& n);
Numerical* ln_light(const Numerical& n);
Numerical* ln(const Numerical& n);
Numerical* log(const Numerical& n, const Numerical& a);
Numerical* exp(const Numerical& n);
Numerical* pow(const Numerical& n, const Numerical& a);
Numerical* cos(const Numerical& n);
Numerical* sin(const Numerical& n);
Numerical* tan(const Numerical& n);
Numerical* sec(const Numerical& n);
Numerical* csc(const Numerical& n);
Numerical* cot(const Numerical& n);
Numerical* arccos(const Numerical& n);
Numerical* arcsin(const Numerical& n);
Numerical* arctan(const Numerical& n);
Numerical* arcsec(const Numerical& n);
Numerical* arccsc(const Numerical& n);
Numerical* arccot(const Numerical& n);
Numerical* cosh(const Numerical& n);
Numerical* sinh(const Numerical& n);
Numerical* tanh(const Numerical& n);
Numerical* sech(const Numerical& n);
Numerical* csch(const Numerical& n);
Numerical* coth(const Numerical& n);
Numerical* arccosh(const Numerical& n);
Numerical* arcsinh(const Numerical& n);
Numerical* arctanh(const Numerical& n);
Numerical* arcsech(const Numerical& n);
Numerical* arccsch(const Numerical& n);
Numerical* arccoth(const Numerical& n);


#endif
