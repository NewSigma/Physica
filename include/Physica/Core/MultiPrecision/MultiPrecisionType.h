/*
 * Copyright 2020-2021 WeiBo He.
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
#pragma once

#include <cstdint>
#include "Physica/SystemBits.h"
/**
 * This header contains some types and definitions used by package MultiPrecision.
 * Here MP stands for MultiPrecision.
 */
#ifdef PHYSICA_64BIT
    #define MPUnitWidth (64U)
    using MPUnit = uint64_t;
    using SignedMPUnit = int64_t;
    #define MPUnitMax INT64_MAX
#endif

#ifdef PHYSICA_32BIT
    #define MPUnitWidth (32U)
    using MPUnit = uint32_t;
    using SignedMPUnit = int32_t;
    #define MPUnitMax INT32_MAX
#endif

#define MPUnitHighestBitMask ((MPUnit)1 << (MPUnitWidth - 1))
#define MPUnitLowerMask (MPUnitMax >> (MPUnitWidth / 2))

namespace Physica::Core {
    enum ScalarOption {
        Float = 0,
        Double = 1,
        MultiPrecision = 2
    };

    /**
     * \class Scalar is a advanced float type that supports multiple precision and error track,
     * which is also compatible with float and double.
     */
    template<ScalarOption option = MultiPrecision, bool errorTrack = true> class Scalar;
}
