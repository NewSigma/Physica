/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_SYSTEMBITS_H
#define PHYSICA_SYSTEMBITS_H

#include <climits>

#if LONG_WIDTH == 64
#define PHYSICA_64BIT
#elif LONG_WIDTH == 32;
#define PHYSICA_32BIT
#else
#error Uncompatible operating system bits.
#endif

#endif //PHYSICA_SYSTEMBITS_H
