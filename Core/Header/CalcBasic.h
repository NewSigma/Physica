/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_CALCBASIC_H
#define PHYSICA_CALCBASIC_H

class Numerical;

Numerical add(const Numerical& n1, const Numerical& n2);
Numerical sub(const Numerical& n1, const Numerical& n2);
Numerical mul(const Numerical& n1, const Numerical& n2);
Numerical div(const Numerical& n1, const Numerical& n2);
bool cutLength(Numerical& n);
void cutZero(Numerical& n);

#endif
