/*
 * Copyright 2019-2024 Weibo He.
 *
 * This file is part of Physica.
 *
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

#include <cmath>
#include <ostream>
#include <random>
#ifdef PHYSICA_CUDA
    #include <cuda/std/limits>
#endif
#include <Physica/Core/Exception/NoImplException.h>
#include <Physica/Core/IO/HDF5/HDF5.h>
#include "MultiPrecisionType.h"
#include "ScalarImpl/ScalarBase.h"

namespace Physica::Core {
    namespace Internal {
        /**
         * This class return a type that can exactly represent the two input scalars.
         */
        template<class AnyScalar1, class AnyScalar2>
        class BinaryScalarOpReturnType {
            static_assert(Core::is_scalar<AnyScalar1>::value, "[Error]: This is not a scalar type");
            static_assert(Core::is_scalar<AnyScalar2>::value, "[Error]: This is not a scalar type");
            constexpr static ScalarOption Option = std::max(Traits<AnyScalar1>::Option, Traits<AnyScalar2>::Option);
            constexpr static bool isComplex = AnyScalar1::isComplex || AnyScalar2::isComplex;
            constexpr static bool isDiffable1 = AnyScalar1::isDifferentiable;
            constexpr static bool isDiffable2 = AnyScalar2::isDifferentiable;
            constexpr static bool isDiffable = isDiffable1 || isDiffable2;

            constexpr static DiffMode Mode = AnyScalar1::isDifferentiable ? Traits<AnyScalar1>::Mode : Traits<AnyScalar2>::Mode;
            constexpr static DiffMode Mode1 = AnyScalar2::isDifferentiable ? Traits<AnyScalar2>::Mode : Traits<AnyScalar1>::Mode;
            static_assert(Mode == Mode1, "[Error]: Operation between differentiable scalars with different mode is not well defined");

            constexpr static int Order1 = Traits<AnyScalar1>::Order;
            constexpr static int Order2 = Traits<AnyScalar2>::Order;
            constexpr static bool UseMixOrder = isDiffable1 && isDiffable2 && (Order1 != Order2);
            constexpr static int Order = UseMixOrder ? std::min(Order1, Order2) : std::max(Order1, Order2);
            static_assert(!(Mode == DiffMode::Reverse && UseMixOrder), "[Error]: Reverse mode does not support mixed order");

            using Type0 = Scalar<Option>;
            using Type1 = typename std::conditional<isComplex, Complex<Type0>, Type0>::type;
            using Type2 = typename std::conditional<isDiffable, Diff<Type1, Mode, Order>, Type1>::type;
        public:
            using Type = Type2;
        };
    }

    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed> class Random;
    template<class ScalarType, class RandomType> class GaussRandomPool;
    template<ScalarOption Option> __host__ __device__ inline Scalar<Option> abs(const Scalar<Option>& s) noexcept;
    template<ScalarOption Option> __host__ __device__ inline Scalar<Option> square(const Scalar<Option>& s) noexcept;
    template<ScalarOption Option> __host__ __device__ inline Scalar<Option> sqrt(const Scalar<Option>& s) noexcept;
    template<ScalarOption Option> __host__ __device__ inline Scalar<Option> ln(const Scalar<Option>& s) noexcept;

    inline int countLeadingZeros(MPUnit n) noexcept;
    inline int countBackZeros(unsigned long n) noexcept;
    inline int countOnes(MPUnit n) noexcept;
    inline double convertDoubleImpl(int length, int power, MPUnit* __restrict byte);
}

#include "Rational.h"
#include "ScalarImpl/ScalarImpl.h"
#ifdef PHYSICA_CUDA
    #include "ScalarImpl/Float16.h"
#endif
#include "ScalarImpl/Float32.h"
#include "ScalarImpl/Float64.h"
#include "ScalarImpl/Math.h"
#include "ScalarImpl/FloatMP.h"
#include "ScalarImpl/SIMD.h"
