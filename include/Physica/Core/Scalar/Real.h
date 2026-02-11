/*
 * Copyright 2019-2026 Weibo He.
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

#include <cstring>
#include <cmath>
#include <bit>
#ifdef PHYSICA_CUDA
    #include <cuda/std/limits>
#endif
#include "Scalar.h"

namespace Physica {
    template<Scalar, class RandomSource> class GaussRandomPool;
    using float16 = Real<Float16>;
    using float32 = Real<Float32>;
    using float64 = Real<Float64>;
    /**
     * We promote to the minimal type if precisions do not match.
     */
    template<FloatPrec Prec1, FloatPrec Prec2>
    [[nodiscard]] __host__ __device__ auto operator+(Real<Prec1> x, Real<Prec2> y) noexcept requires(Prec1 != Prec2) {
        using ResultType = Internal::BinaryScalarOpRtnTy<Real<Prec1>, Real<Prec2>>::Type;
        return ResultType(x) + ResultType(y);
    }

    template<FloatPrec Prec1, FloatPrec Prec2>
    [[nodiscard]] __host__ __device__ auto operator-(Real<Prec1> x, Real<Prec2> y) noexcept requires(Prec1 != Prec2) {
        using ResultType = Internal::BinaryScalarOpRtnTy<Real<Prec1>, Real<Prec2>>::Type;
        return ResultType(x) - ResultType(y);
    }

    template<FloatPrec Prec1, FloatPrec Prec2>
    [[nodiscard]] __host__ __device__ auto operator*(Real<Prec1> x, Real<Prec2> y) noexcept requires(Prec1 != Prec2) {
        using ResultType = Internal::BinaryScalarOpRtnTy<Real<Prec1>, Real<Prec2>>::Type;
        return ResultType(x) * ResultType(y);
    }

    template<FloatPrec Prec1, FloatPrec Prec2>
    [[nodiscard]] __host__ __device__ auto operator/(Real<Prec1> x, Real<Prec2> y) noexcept requires(Prec1 != Prec2) {
        using ResultType = Internal::BinaryScalarOpRtnTy<Real<Prec1>, Real<Prec2>>::Type;
        return ResultType(x) / ResultType(y);
    }

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

    template<FloatPrec Prec>
    auto calcULPDiff(Real<Prec> x1, Real<Prec> x2) noexcept {
        using PossibleRtnTy = std::tuple<uint16_t, uint32_t, uint64_t>;
        using Int = std::remove_reference<decltype(std::get<Prec>(std::declval<PossibleRtnTy>()))>::type;
        constexpr Int SignMask = Int(1) << (sizeof(Int) * CHAR_BIT - 1);
        constexpr Int AntiSignMask = ~SignMask;
        static_assert(sizeof(Int) == sizeof(Real<Prec>));

        assert(x1.isFinite() && x2.isFinite());
        Int i1, i2;
        memcpy(&i1, &x1, sizeof(Int));
        memcpy(&i2, &x2, sizeof(Int));
        bool signdiff = (i1 ^ i2) & SignMask;
        if (signdiff)
            return (i1 & AntiSignMask) + (i2 & AntiSignMask);
        return i1 > i2 ? (i1 - i2) : (i2 - i1);
    }
}

#include "Rational.h"
#include "Physica/Core/IO/HDF5/HDF5.h"
#include "Physica/Core/Math/Random/Random.h"
#include "ScalarImpl/ScalarBase.h"
#include "RealImpl/Float32.h"
#include "RealImpl/Float64.h"
#include "RealImpl/FloatMP.h"
#include "RealImpl/MathConst.h"
#include "RealImpl/Math.h"
#ifdef PHYSICA_CUDA 
    #include "RealImpl/Float16.h"
    #include "RealImpl/MathFP16.h"
#endif
#include "RealImpl/MathFPMP.h"
#include "RealImpl/SIMD.h"
