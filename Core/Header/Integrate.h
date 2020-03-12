#ifndef PHYSICA_INTEGRATE_H
#define PHYSICA_INTEGRATE_H

#include "AbstractNum.h"

enum IntegrateMethod {
    Rectangular,
    Ladder
};

AbstractNum* Integrate(AbstractNum* func(AbstractNum*), AbstractNum* x0, AbstractNum* x1, IntegrateMethod method);
AbstractNum* rectangular(AbstractNum* func(AbstractNum*), AbstractNum* x0, AbstractNum* x1);
AbstractNum* ladder(AbstractNum* func(AbstractNum*), AbstractNum* x0, AbstractNum* x1);

#endif
