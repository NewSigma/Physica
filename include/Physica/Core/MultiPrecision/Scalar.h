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
#include <concepts>
#include "Physica/Macro.h"

namespace Physica::Core {
    template<class Derived> class ScalarBase;

    template<class T>
    concept Scalar = std::derived_from<T, ScalarBase<T>>;

    enum ScalarOption {
        Float16 = 0,
        Float32 = 1,
        Float64 = 2,
        FloatMP = 3,

        Half = Float16,
        Float = Float32,
        Double = Float64
    };

    enum class DiffMode {
        Forward,
        Reverse
    };

    // MP = MultiPrecision
    using MPUnit = std::conditional<PhysicaWordSize == 64, uint64_t, uint32_t>::type;
    using SignedMPUnit = std::conditional<PhysicaWordSize == 64, int64_t, int32_t>::type;
    constexpr static unsigned int MPUnitWidth = PhysicaWordSize;
    constexpr static MPUnit MPUnitMax = std::numeric_limits<MPUnit>::max();
    constexpr static MPUnit MPUnitHighestBitMask = static_cast<MPUnit>(1) << (MPUnitWidth - 1);
    constexpr static MPUnit MPUnitLowerMask = MPUnitMax >> (MPUnitWidth / 2);
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

    template<class ScalarType, DiffMode Mode, int Order = 1>
    class Diff;

    namespace Internal {
        /**
         * \class BinaryScalarOpRtnTy returns a type that can exactly represent the two input scalars.
         */
        template<Scalar T1, Scalar T2>
        class BinaryScalarOpRtnTy {
            constexpr static ScalarOption Option = std::max(T1::Option, T2::Option);
            constexpr static bool isComplex = T1::isComplex || T2::isComplex;
            constexpr static bool isDiffable1 = T1::isDifferentiable;
            constexpr static bool isDiffable2 = T2::isDifferentiable;
            constexpr static bool isDiffable = isDiffable1 || isDiffable2;

            constexpr static DiffMode Mode = T1::isDifferentiable ? T1::Mode : T2::Mode;
            constexpr static DiffMode Mode1 = T2::isDifferentiable ? T2::Mode : T1::Mode;
            static_assert(Mode == Mode1, "[Error]: Operation between differentiable scalars with different mode is not well defined");

            constexpr static int Order1 = T1::Order;
            constexpr static int Order2 = T2::Order;
            constexpr static bool UseMixOrder = isDiffable1 && isDiffable2 && (Order1 != Order2);
            constexpr static int Order = UseMixOrder ? std::min(Order1, Order2) : std::max(Order1, Order2);
            static_assert(!(Mode == DiffMode::Reverse && UseMixOrder), "[Error]: Reverse mode does not support mixed order");

            using Type0 = Real<Option>;
            using Type1 = std::conditional<isComplex, Complex<Type0>, Type0>::type;
            using Type2 = std::conditional<isDiffable, Diff<Type1, Mode, Order>, Type1>::type;
        public:
            using Type = Type2;
        };
    }
}
