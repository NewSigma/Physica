/*
 * Copyright 2020-2024 WeiBo He.
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
#include "Physica/Macro.h"
/**
 * This header contains some types and definitions used by package MultiPrecision.
 * Here MP stands for MultiPrecision.
 */
namespace Physica::Core {
    using MPUnit = typename std::conditional<PhysicaWordSize == 64, uint64_t, uint32_t>::type;
    using SignedMPUnit = typename std::conditional<PhysicaWordSize == 64, int64_t, int32_t>::type;
    constexpr static unsigned int MPUnitWidth = PhysicaWordSize;
    constexpr static MPUnit MPUnitMax = std::numeric_limits<MPUnit>::max();
    constexpr static MPUnit MPUnitHighestBitMask = static_cast<MPUnit>(1) << (MPUnitWidth - 1);
    constexpr static MPUnit MPUnitLowerMask = MPUnitMax >> (MPUnitWidth / 2);

    enum ScalarOption {
        Float = 0,
        Double = 1,
        MultiPrecision = 2
    };

    /**
     * \class Scalar is a advanced float type that supports multiple precision
     */
    template<ScalarOption Option = Double> class Scalar;
    template<class AnyScalar> class ComplexScalar;
}
