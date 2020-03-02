#ifndef PHYSICA_INTEGRATE_H
#define PHYSICA_INTEGRATE_H

#include "Number.h"

enum IntegrateMethod {
    Rectangular,
    Ladder
};

Number* Integrate(Number* func(Number*), Number* x0, Number* x1, IntegrateMethod method);
Number* rectangular(Number* func(Number*), Number* x0, Number* x1);
Number* ladder(Number* func(Number*), Number* x0, Number* x1);

#endif
