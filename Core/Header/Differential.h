#ifndef PHYSICA_DIFFERENTIAL_H
#define PHYSICA_DIFFERENTIAL_H

#include "RealNumber.h"

enum DifferentialMethod{
    Identify
};

RealNumber* D_double_point(RealNumber* func(const RealNumber&), const RealNumber& x0);
RealNumber* D_right(RealNumber* func(const RealNumber&), const RealNumber& x0);
RealNumber* D_left(RealNumber* func(const RealNumber&), const RealNumber& x0);

#endif
