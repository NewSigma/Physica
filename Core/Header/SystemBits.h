/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_SYSTEMBITS_H
#define PHYSICA_SYSTEMBITS_H

#include <climits>

#ifdef __WORDSIZE
    #if __WORDSIZE == 64
        #define PHYSICA_64BIT
    #elif __WORDSIZE == 32;
        #define PHYSICA_32BIT
    #else
        #error Uncompatible operating system bits.
#endif
#else
    #error __WORDSIZE not defined.
#endif

#if __WORDSIZE == INT_WIDTH
typedef unsigned int NumericalUnit;
typedef int SignedNumericalUnit;
const NumericalUnit NumericalUnitMax = UINT_MAX;
const size_t NumericalUnitWidth = INT_WIDTH;
#elif __WORDSIZE == LONG_WIDTH
typedef unsigned long NumericalUnit;
typedef long SignedNumericalUnit;
const NumericalUnit NumericalUnitMax = ULONG_MAX;
const size_t NumericalUnitWidth = LONG_WIDTH;
#elif __WORDSIZE == LLONG_WIDTH
typedef unsigned long long NumericalUnit;
typedef long long SignedNumericalUnit;
const NumericalUnit NumericalUnitMax = ULLLONG_MAX;
const size_t NumericalUnitWidth = LLONG_WIDTH;
#else
#error No marching NumericalUnit.
#endif

#endif
