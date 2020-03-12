#ifndef PHYSICA_DIFFERENTIAL_H
#define PHYSICA_DIFFERENTIAL_H

#include "AbstractNum.h"

AbstractNum* D_double_point(AbstractNum* func(const AbstractNum&), const AbstractNum& x0);
AbstractNum* D_right(AbstractNum* func(const AbstractNum&), const AbstractNum& x0);
AbstractNum* D_left(AbstractNum* func(const AbstractNum&), const AbstractNum& x0);

#endif
