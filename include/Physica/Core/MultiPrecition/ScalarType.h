/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_SCALARTYPE_H
#define PHYSICA_SCALARTYPE_H

#include "Physica/SystemBits.h"
/*!
 * This header contains some types and definitions used by Scalar.
 */
#ifdef PHYSICA_64BIT
    #define ScalarUnitWidth (64U)
#endif

#ifdef PHYSICA_32BIT
    #define ScalarUnitWidth (32U)
#endif

#if PhysicaWordSize == INT_WIDTH
    typedef unsigned int ScalarUnit;
    typedef int SignedScalarUnit;
    #define ScalarUnitMax UINT_MAX
#elif PhysicaWordSize == LONG_WIDTH
    typedef unsigned long ScalarUnit;
    typedef long SignedScalarUnit;
    #define ScalarUnitMax ULONG_MAX
#elif PhysicaWordSize == LLONG_WIDTH
    typedef unsigned long long ScalarUnit;
    typedef long long SignedScalarUnit;
    #define ScalarUnitMax ULLLONG_MAX
#else
    #error No marching ScalarUnit.
#endif

#define ScalarUnitHighestBitMask ((ScalarUnit)1 << (ScalarUnitWidth - 1))
#define ScalarUnitLowerMask (ScalarUnitMax >> (ScalarUnitWidth / 2))

//Use of double_extract may cause several warnings in valgrind.
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