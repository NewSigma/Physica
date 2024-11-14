/*
 * Copyright 2020-2024 Weibo He.
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
#include <limits>
#include <type_traits>
#include <Physica/Macro.h>

namespace Physica::Core {
    template<class T> class ScalarBase;

    template<class T>
    struct is_scalar : public std::is_base_of<ScalarBase<T>, T> {};

    template<class T>
    struct is_scalar<ScalarBase<T>> {
        constexpr static bool value = true;
    };

    template<class T>
    concept Scalar = is_scalar<T>::value;
    // MP = MultiPrecision
    using MPUnit = typename std::conditional<PhysicaWordSize == 64, uint64_t, uint32_t>::type;
    using SignedMPUnit = typename std::conditional<PhysicaWordSize == 64, int64_t, int32_t>::type;
    constexpr static unsigned int MPUnitWidth = PhysicaWordSize;
    constexpr static MPUnit MPUnitMax = std::numeric_limits<MPUnit>::max();
    constexpr static MPUnit MPUnitHighestBitMask = static_cast<MPUnit>(1) << (MPUnitWidth - 1);
    constexpr static MPUnit MPUnitLowerMask = MPUnitMax >> (MPUnitWidth / 2);

    enum ScalarOption {
        Float16 = 0,
        Float32 = 1,
        Float64 = 2,
        FloatMP = 3,

        Half = Float16,
        Float = Float32,
        Double = Float64
    };
    /**
     * \class Real is a advanced float type that supports multiple precision
     */
    template<ScalarOption Option = Float64> class Real;
    using float16 = Real<Float16>;
    using float32 = Real<Float32>;
    using float64 = Real<Float64>;

    template<class T> class Complex;
    using cfloat16 = Complex<float16>;
    using cfloat32 = Complex<float32>;
    using cfloat64 = Complex<float64>;
}
