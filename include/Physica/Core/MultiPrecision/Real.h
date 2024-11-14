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
#include "Scalar.h"
#include "RealImpl/ScalarBase.h"

namespace Physica::Core {
    namespace Internal {
        /**
         * This class return a type that can exactly represent the two input scalars.
         */
        template<Scalar T1, Scalar T2>
        class BinaryScalarOpReturnType {
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
            using Type1 = typename std::conditional<isComplex, Complex<Type0>, Type0>::type;
            using Type2 = typename std::conditional<isDiffable, Diff<Type1, Mode, Order>, Type1>::type;
        public:
            using Type = Type2;
        };
    }

    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed> class Random;
    template<class ScalarType, class RandomType> class GaussRandomPool;
    template<ScalarOption Option> __host__ __device__ inline Real<Option> abs(const Real<Option>& s) noexcept;
    template<ScalarOption Option> __host__ __device__ inline Real<Option> square(const Real<Option>& s) noexcept;
    template<ScalarOption Option> __host__ __device__ inline Real<Option> sqrt(const Real<Option>& s) noexcept;
    template<ScalarOption Option> __host__ __device__ inline Real<Option> ln(const Real<Option>& s) noexcept;

    inline int countLeadingZeros(MPUnit n) noexcept;
    inline int countBackZeros(unsigned long n) noexcept;
    inline int countOnes(MPUnit n) noexcept;
    inline double convertDoubleImpl(int length, int power, MPUnit* __restrict byte);
}

#include "Rational.h"
#include "RealImpl/RealImpl.h"
#ifdef PHYSICA_CUDA
    #include "RealImpl/Float16.h"
#endif
#include "RealImpl/Float32.h"
#include "RealImpl/Float64.h"
#include "RealImpl/Math.h"
#include "RealImpl/FloatMP.h"
#include "RealImpl/SIMD.h"
