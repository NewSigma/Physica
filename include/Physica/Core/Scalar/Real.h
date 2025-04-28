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
#include <bit>
#ifdef PHYSICA_CUDA
    #include <cuda/std/limits>
#endif
#include "Scalar.h"

namespace Physica {
    template<Scalar, class RandomSource> class GaussRandomPool;

    union float_extract {
    private:
        constexpr static bool Flag = std::endian::native == std::endian::little;
    public:
        float value;
        struct {
            unsigned int low : Flag ? 16 : 1;
            unsigned int high : Flag ? 7 : 8;
            unsigned int exp : Flag ? 8 : 7;
            unsigned int sign : Flag ? 1 : 16;
        };
    };

    union double_extract {
    private:
        constexpr static bool Flag = std::endian::native == std::endian::little;
    public:
        double value;
        struct {
            unsigned int low : Flag ? 32 : 1;
            unsigned int high : Flag ? 20 : 11;
            unsigned int exp : Flag ? 11 : 20;
            unsigned int sign : Flag ? 1 : 32;
        };
    };

    template<FloatPrec Prec1, FloatPrec Prec2>
    __host__ __device__ auto operator+(Real<Prec1> x, Real<Prec2> y) requires(Prec1 != Prec2) {
        using ResultType = Internal::BinaryScalarOpRtnTy<Real<Prec1>, Real<Prec2>>::Type;
        return ResultType(x) + ResultType(y);
    }

    template<FloatPrec Prec1, FloatPrec Prec2>
    __host__ __device__ auto operator-(Real<Prec1> x, Real<Prec2> y) requires(Prec1 != Prec2) {
        using ResultType = Internal::BinaryScalarOpRtnTy<Real<Prec1>, Real<Prec2>>::Type;
        return ResultType(x) - ResultType(y);
    }

    template<FloatPrec Prec1, FloatPrec Prec2>
    __host__ __device__ auto operator*(Real<Prec1> x, Real<Prec2> y) requires(Prec1 != Prec2) {
        using ResultType = Internal::BinaryScalarOpRtnTy<Real<Prec1>, Real<Prec2>>::Type;
        return ResultType(x) * ResultType(y);
    }

    template<FloatPrec Prec1, FloatPrec Prec2>
    __host__ __device__ auto operator/(Real<Prec1> x, Real<Prec2> y) requires(Prec1 != Prec2) {
        using ResultType = Internal::BinaryScalarOpRtnTy<Real<Prec1>, Real<Prec2>>::Type;
        return ResultType(x) / ResultType(y);
    }

    inline double convertDoubleImpl(int length, int power, MPUnit* __restrict byte);

    using float16 = Real<Float16>;
    using float32 = Real<Float32>;
    using float64 = Real<Float64>;
}

#include "Rational.h"
#include "RealImpl/RealImpl.h"

#include "Physica/Core/IO/HDF5/HDF5.h"
#include "Physica/Core/Math/Random/Random.h"
#include "ScalarImpl/ScalarBase.h"
#include "RealImpl/Float32.h"
#include "RealImpl/Float64.h"
#ifdef PHYSICA_CUDA
    #include "RealImpl/Float16.h"
#endif
#include "RealImpl/FloatMP.h"
#include "RealImpl/MathConst.h"
#include "RealImpl/Math.h"
#ifdef PHYSICA_CUDA 
    #include "RealImpl/MathFP16.h"
#endif
#include "RealImpl/MathFPMP.h"
#include "RealImpl/SIMD.h"
