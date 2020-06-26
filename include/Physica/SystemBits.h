/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_SYSTEMBITS_H
#define PHYSICA_SYSTEMBITS_H

#include <climits>
#include "Physica/Config.h"

#if PhysicaWordSize == 64
    #define PHYSICA_64BIT
#elif PhysicaWordSize == 32;
    #define PHYSICA_32BIT
    #define ScalarUnitWidth (32U)
#endif

#endif
