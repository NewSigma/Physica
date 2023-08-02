/*
 * Copyright 2020-2023 WeiBo He.
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

#include <iomanip>
#include <random>
#include "ScalarArithmetic.h"
/**
 * This file is part of implementations of \Scalar.
 * Do not include this header file, include \file Scalar.h instead.
 */
namespace Physica::Core {
    //////////////////////////////////////////////Global//////////////////////////////////////////////
    template<class ScalarType>
    ScalarType relativeError(const ScalarType& scalar1, const ScalarType& scalar2) {
        static_assert(Core::is_scalar<ScalarType>::value && !ScalarType::isComplex && !ScalarType::isDifferentiable, "[Error]: Invalid template param");
        const auto& s1 = scalar1.getDerived();
        const auto& s2 = scalar2.getDerived();
        const ScalarType min = std::numeric_limits<ScalarType>::min();
        const bool useAbsCompare = (abs(s1) < min) || (abs(s2) < min);
        const ScalarType delta = s1 - s2;
        const ScalarType error = useAbsCompare ? abs(delta) : abs(delta / (s1 + s2) * ScalarType(2));
        return error;
    }

    template<class ScalarType>
    bool scalarNear(const ScalarBase<ScalarType>& s1,
                    const ScalarBase<ScalarType>& s2,
                    double precision) {
        assert(precision > 0);
        constexpr bool isDifferentiable = ScalarType::isDifferentiable;
        if constexpr (ScalarType::isComplex) {
            using PlainRealType = typename ScalarType::PlainScalar::RealType;
            const ScalarType diff = s1.getDerived() - s2.getDerived();
            const bool isValueNear = scalarNear(abs(diff.getValue()), PlainRealType(0), precision);
            if constexpr (isDifferentiable)
                return isValueNear && scalarNear(abs(diff.getTangent()), PlainRealType(0), precision);
            else
                return isValueNear;
        }
        else {
            using PlainScalar = typename ScalarType::PlainScalar;
            const bool isValueNear = relativeError(s1.getValue().getReal(), s2.getValue().getReal()) < PlainScalar(precision);
            if constexpr (ScalarType::isDifferentiable)
                return isValueNear && relativeError(s1.getTangent().getReal(), s2.getTangent().getReal()) < PlainScalar(precision);
            else
                return isValueNear;
        }
    }

    template<ScalarOption option>
    std::ostream& operator<<(std::ostream& os, const Scalar<option>& s) {
        return os << std::setprecision(10) //10 is the max precision of double.
                  << double(s)
                  << std::setprecision(6); //6 is the default precision.
    }

    template<ScalarOption option>
    inline Scalar<option> operator+(const Scalar<option>& s) {
        return Scalar<option>(s);
    }

    template<class ScalarType1, class ScalarType2>
    __host__ __device__ inline void operator+=(ScalarBase<ScalarType1>& s1, const ScalarBase<ScalarType2>& s2) {
        s1.getDerived() = s1.getDerived() + s2.getDerived();
    }

    template<class ScalarType1, class ScalarType2>
    __host__ __device__ inline void operator-=(ScalarBase<ScalarType1>& s1, const ScalarBase<ScalarType2>& s2) {
        s1.getDerived() = s1.getDerived() - s2.getDerived();
    }

    template<class ScalarType1, class ScalarType2>
    __host__ __device__ inline void operator*=(ScalarBase<ScalarType1>& s1, const ScalarBase<ScalarType2>& s2) {
        s1.getDerived() = s1.getDerived() * s2.getDerived();
    }

    template<class ScalarType1, class ScalarType2>
    __host__ __device__ inline void operator/=(ScalarBase<ScalarType1>& s1, const ScalarBase<ScalarType2>& s2) {
        s1.getDerived() = s1.getDerived() / s2.getDerived();
    }

    template<ScalarOption option>
    inline void operator^=(Scalar<option>& s1, const Scalar<option>& s2) { s1 = s1 ^ s2; }

    template<ScalarOption option>
    inline void operator<<=(Scalar<option>& s, int bits) { s = s << bits; }

    template<ScalarOption option>
    inline void operator>>=(Scalar<option>& s, int bits) { s = s >> bits; }

    template<ScalarOption option>
    __host__ __device__ inline bool operator>=(const Scalar<option>& s1, const Scalar<option>& s2) {
        return !(s1 < s2);
    }

    template<ScalarOption option>
    __host__ __device__ inline bool operator<=(const Scalar<option>& s1, const Scalar<option>& s2) {
        return !(s1 > s2);
    }

    template<ScalarOption option>
    __host__ __device__ inline bool operator!= (const Scalar<option>& s1, const Scalar<option>& s2) {
        return !(s1 == s2);
    }

    template<ScalarOption option>
    inline void swap(Scalar<option>& s1, Scalar<option>& s2) noexcept {
        s1.swap(s2);
    }
    ///////////////////////////////////////////MultiPrecision/////////////////////////////////////////
    /**
     * Returns true if s1 and s2 has the same sign. Both s1 and s2 do not equal to zero.
     * This function provide a quick sign check compare to using isPositive() and isNegative().
     */
    inline bool Scalar<MultiPrecision>::matchSign(const Scalar<MultiPrecision>& s1, const Scalar<MultiPrecision>& s2) {
        assert(!s1.isZero() && !s2.isZero());
        return (s1.length ^ s2.length) >= 0; //NOLINT Bitwise operator between two signed integer is intended.
    }
    /**
     * Cut zeros from the beginning.
     */
    inline void Scalar<MultiPrecision>::cutZero() {
        const int size = getSize();
        int id = size - 1;
        while(byte[id] == 0 && id > 0)
            --id;
        ++id;

        if(id != size) {
            int shorten = size - id;
            byte = reinterpret_cast<MPUnit*>(realloc(byte, id * sizeof(MPUnit)));
            length = length > 0 ? id : -id;
            auto temp = power;
            power = byte[id - 1] != 0 ? (temp - shorten) : 0;
        }
    }

    inline Scalar<MultiPrecision>& operator++(Scalar<MultiPrecision>& s) {
        s += BasicConst::getInstance()._1;
        return s;
    }

    inline Scalar<MultiPrecision>& operator--(Scalar<MultiPrecision>& s) {
        s -= BasicConst::getInstance()._1;
        return s;
    }

    inline Scalar<MultiPrecision> operator++(Scalar<MultiPrecision>& s, int) {
        Scalar<MultiPrecision> temp(s);
        s += BasicConst::getInstance()._1;
        return temp;
    }

    inline Scalar<MultiPrecision> operator--(Scalar<MultiPrecision>& s, int) {
        Scalar<MultiPrecision> temp(s);
        s -= BasicConst::getInstance()._1;
        return temp;
    }
    /////////////////////////////////////////////Float////////////////////////////////////////////////
    __host__ __device__ inline Scalar<Float>::Scalar(const Scalar<Double>& s) : f(float(s)) {}

    inline Scalar<Float>& operator++(Scalar<Float>& s) {
        s += Scalar<Float>(1.0F);
        return s;
    }

    inline Scalar<Float>& operator--(Scalar<Float>& s) {
        s -= Scalar<Float>(1.0F);
        return s;
    }

    inline Scalar<Float> operator++(Scalar<Float>& s, int) {
        Scalar<Float> temp(s);
        s += Scalar<Float>(1.0F);
        return temp;
    }

    inline Scalar<Float> operator--(Scalar<Float>& s, int) {
        Scalar<Float> temp(s);
        s -= Scalar<Float>(1.0F);
        return temp;
    }

    inline bool absCompare(const Scalar<Float>& s1, const Scalar<Float>& s2) {
        return fabsf(s1.getTrivial()) >= fabsf(s2.getTrivial());
    }

    __host__ __device__ inline bool operator> (const Scalar<Float>& s1, const Scalar<Float>& s2) {
        return s1.getTrivial() > s2.getTrivial();
    }

    __host__ __device__ inline bool operator< (const Scalar<Float>& s1, const Scalar<Float>& s2) {
        return s1.getTrivial() < s2.getTrivial();
    }

    __host__ __device__ inline bool operator== (const Scalar<Float>& s1, const Scalar<Float>& s2) {
        return s1.getTrivial() == s2.getTrivial();
    }

    template<class RandomGenerator>
    Scalar<Float> Scalar<Float>::random_uniform(RandomGenerator& gen) {
        std::uniform_real_distribution<float> dist{};
        return Scalar(dist(gen));
    }

    template<class RandomGenerator>
    Scalar<Float> Scalar<Float>::random_normal(RandomGenerator& gen) {
        std::normal_distribution<float> dist{};
        return Scalar(dist(gen));
    }
    /////////////////////////////////////////////Double////////////////////////////////////////////////
    __host__ __device__ inline Scalar<Double>::Scalar(const Scalar<Float>& s) : d(double(s)) {}

    inline Scalar<Double>& operator++(Scalar<Double>& s) {
        s += Scalar<Double>(1.0);
        return s;
    }

    inline Scalar<Double>& operator--(Scalar<Double>& s) {
        s -= Scalar<Double>(1.0);
        return s;
    }

    inline Scalar<Double> operator++(Scalar<Double>& s, int) {
        Scalar<Double> temp(s);
        s += Scalar<Double>(1.0);
        return temp;
    }

    inline Scalar<Double> operator--(Scalar<Double>& s, int) {
        Scalar<Double> temp(s);
        s -= Scalar<Double>(1.0);
        return temp;
    }

    inline bool absCompare(const Scalar<Double>& s1, const Scalar<Double>& s2) {
        return fabs(s1.getTrivial()) >= fabs(s2.getTrivial());
    }

    __host__ __device__ inline bool operator> (const Scalar<Double>& s1, const Scalar<Double>& s2) {
        return s1.getTrivial() > s2.getTrivial();
    }

    __host__ __device__ inline bool operator< (const Scalar<Double>& s1, const Scalar<Double>& s2) {
        return s1.getTrivial() < s2.getTrivial();
    }

    __host__ __device__ inline bool operator== (const Scalar<Double>& s1, const Scalar<Double>& s2) {
        return s1.getTrivial() == s2.getTrivial();
    }

    template<class RandomGenerator>
    Scalar<Double> Scalar<Double>::random_uniform(RandomGenerator& gen) {
        std::uniform_real_distribution<double> dist{};
        return Scalar(dist(gen));
    }

    template<class RandomGenerator>
    Scalar<Double> Scalar<Double>::random_normal(RandomGenerator& gen) {
        std::normal_distribution<double> dist{};
        return Scalar(dist(gen));
    }
}

#include "ScalarArithmetic.h"
#include "ElementaryFunction.h"
#include "Operation/Pow.h"

namespace Physica::Core {
    template<ScalarOption option>
    inline Scalar<option> operator^(const Scalar<option>& s1, const Scalar<option>& s2) {
        return pow(s1, s2);
    }

    template<>
    inline Scalar<MultiPrecision> operator^(const Scalar<MultiPrecision>& s1, const Scalar<MultiPrecision>& s2) {
        return s1.isInteger() ? powScalar(s1, s2) : exp(ln(s1) * s2);
    }
}
