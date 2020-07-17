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

/*!
 * Improve: Platform dependent, may be wrong on some spatial platforms, add a test before compile.
 * Note: Use of double_extract may cause several warnings in valgrind.
 */
#if PhysicaEndian == PhysicaBigEndian
    union float_extract {
        float value;
        struct {
            unsigned int sign : 1;
            unsigned int exp : 8;
            unsigned int high : 7;
            unsigned int low : 16;
        };
    };

    union double_extract {
        double value;
        struct {
            unsigned int sign : 1;
            unsigned int exp : 11;
            unsigned int high : 20;
            unsigned int low : 32;
        };
    };
#elif PhysicaEndian == PhysicaLittleEndian
    union float_extract {
        float value;
        struct {
            unsigned int low : 16;
            unsigned int high : 7;
            unsigned int exp : 8;
            unsigned int sign : 1;
        };
    };

    union double_extract {
        double value;
        struct {
            unsigned int low : 32;
            unsigned int high : 20;
            unsigned int exp : 11;
            unsigned int sign : 1;
        };
    };
#endif

#endif
