#ifndef PHYSICA_DIFFERENTIAL_H
#define PHYSICA_DIFFERENTIAL_H

#include "Numerical.h"

Numerical* D_double_point(Numerical* func(const Numerical&), const Numerical& x0);
Numerical* D_right(Numerical* func(const Numerical&), const Numerical& x0);
Numerical* D_left(Numerical* func(const Numerical&), const Numerical& x0);

#endif
