#ifndef PHYSICA_BASICCALCULATES_H
#define PHYSICA_BASICCALCULATES_H

#include "RealNumber.h"

RealNumber* reciprocal(const RealNumber* n);
RealNumber* sqrt(const RealNumber* n, bool checkAcc = true);
RealNumber* factorial(const RealNumber* n);
RealNumber* cos(const RealNumber* n);
RealNumber* sin(const RealNumber* n);
RealNumber* tan(const RealNumber* n);
RealNumber* arccos(const RealNumber* n);
RealNumber* arcsin(const RealNumber* n);
RealNumber* ln(const RealNumber* n);
RealNumber* log(const RealNumber* n, const RealNumber* a);

#endif
