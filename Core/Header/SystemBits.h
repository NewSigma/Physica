/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_SYSTEMBITS_H
#define PHYSICA_SYSTEMBITS_H

#include <climits>
#include "PhysicaConfig.h"

#if PhysicaWordSize == 64
#define PHYSICA_64BIT
#define NumericalUnitWidth (64U)
#elif PhysicaWordSize == 32;
#define PHYSICA_32BIT
#define NumericalUnitWidth (32U)
#endif

#if PhysicaWordSize == INT_WIDTH
typedef unsigned int NumericalUnit;
typedef int SignedNumericalUnit;
#define NumericalUnitMax UINT_MAX
#elif PhysicaWordSize == LONG_WIDTH
typedef unsigned long NumericalUnit;
typedef long SignedNumericalUnit;
#define NumericalUnitMax ULONG_MAX
#elif PhysicaWordSize == LLONG_WIDTH
typedef unsigned long long NumericalUnit;
typedef long long SignedNumericalUnit;
#define NumericalUnitMax ULLLONG_MAX
#else
#error No marching NumericalUnit.
#endif

#define numericalUnitHighestBitMask ((NumericalUnit)1 << (NumericalUnitWidth - 1))
#define numericalUnitLowMask (NumericalUnitMax >> (NumericalUnitWidth / 2))

#if PhysicaEndian == PhysicaBigEndian
union double_extract {
            double value;
            struct {
                unsigned int sign : 1;
                unsigned int exp : 11;
                unsigned int high : 20;
                unsigned int low : 32;
            } structure;
        };
#elif PhysicaEndian == PhysicaLittleEndian
union double_extract {
    double value;
    struct {
        unsigned int low : 32;
        unsigned int high : 20;
        unsigned int exp : 11;
        unsigned int sign : 1;
    } structure;
};
#endif

#endif
