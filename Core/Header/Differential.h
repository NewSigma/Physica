#ifndef PHYSICA_DIFFERENTIAL_H
#define PHYSICA_DIFFERENTIAL_H

#include "RealNumber.h"

RealNumber* derivative_right(RealNumber* func(const RealNumber&), const RealNumber& x0);
RealNumber* derivative_left(RealNumber* func(const RealNumber&), const RealNumber& x0);

#endif
