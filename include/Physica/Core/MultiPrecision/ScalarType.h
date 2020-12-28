/*
 * Copyright 2020 WeiBo He.
 *
 * This file is part of Physica.

 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
 */
#ifndef PHYSICA_SCALARTYPE_H
#define PHYSICA_SCALARTYPE_H

#include "Physica/SystemBits.h"
/*!
 * This header contains some types and definitions used by Scalar.
 */
namespace Physica::Core {
    enum ScalarType {
        Float,
        Double,
        MultiPrecision
    };
}

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

#endif