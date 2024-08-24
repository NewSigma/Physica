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
#include "MultiPrecisionType.h"
#include "ScalarImpl/ScalarBase.h"
#include "Physica/Core/Exception/NoImplException.h"
#include "Physica/Core/IO/HDF5/HDF5.h"

namespace Physica::Core {
    namespace Internal {
        /**
         * This class return a type that can exactly represent the two input scalars.
         */
        template<class AnyScalar1, class AnyScalar2>
        class BinaryScalarOpReturnType {
            static_assert(Core::is_scalar<AnyScalar1>::value, "[Error]: This is not a scalar type");
            static_assert(Core::is_scalar<AnyScalar2>::value, "[Error]: This is not a scalar type");
            static_assert(!AnyScalar1::isComplex && !AnyScalar2::isComplex, "[Error]: This class applies to real scalar only");
            static_assert(!AnyScalar1::isDifferentiable && !AnyScalar2::isDifferentiable, "[Error]: This class applies to plain scalar only");
            static constexpr ScalarOption Option =
                    Traits<AnyScalar1>::Option > Traits<AnyScalar2>::Option ? Traits<AnyScalar1>::Option : Traits<AnyScalar2>::Option;
        public:
            using Type = Scalar<Option>;
        };
    }

    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed> class RandomPool;
    template<class ScalarType, class RandomPoolType> class GaussRandomPool;
    template<ScalarOption Option> __host__ __device__ inline Scalar<Option> abs(const Scalar<Option>& s) noexcept;
    template<ScalarOption Option> __host__ __device__ inline Scalar<Option> square(const Scalar<Option>& s) noexcept;
    template<ScalarOption Option> __host__ __device__ inline Scalar<Option> sqrt(const Scalar<Option>& s) noexcept;
    template<ScalarOption Option> __host__ __device__ inline Scalar<Option> ln(const Scalar<Option>& s) noexcept;

    using float16 = Scalar<Float16>;
    using float32 = Scalar<Float32>;
    using float64 = Scalar<Float64>;
}

#include "Rational.h"
#ifdef PHYSICA_CUDA
    #include "ScalarImpl/Float16.h"
#endif
#include "ScalarImpl/Float32.h"
#include "ScalarImpl/Float64.h"
#include "ScalarImpl/MultiPrecision.h"

#include "ScalarImpl/ElementaryFunction.h"
#include "SIMD.h"
