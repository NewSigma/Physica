#ifndef PHYSICA_SOLVE_H
#define PHYSICA_SOLVE_H

#include "RealNumber.h"
RealNumber* bisectionMethod(RealNumber* func(const RealNumber*), const RealNumber* n, const RealNumber* x1, const RealNumber* x2);
RealNumber* bisectionMethod(RealNumber* func(const RealNumber*), const RealNumber* n, const RealNumber* x_left, const RealNumber* x_right, const RealNumber* y_left, const RealNumber* y2);

#endif
