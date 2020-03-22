#ifndef PHYSICA_INTEGRATE_H
#define PHYSICA_INTEGRATE_H

#include "Numerical.h"

Numerical* rectangular(Numerical* func(Numerical*), Numerical* x0, Numerical* x1);
Numerical* ladder(Numerical* func(Numerical*), Numerical* x0, Numerical* x1);
Numerical* Simpson(Numerical* func(Numerical*), Numerical* x0, Numerical* x1);

#endif
