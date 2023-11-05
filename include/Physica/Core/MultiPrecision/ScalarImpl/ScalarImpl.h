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

    template<class RandomGenerator>
    inline Scalar<Float> Scalar<Float>::random_uniform(RandomGenerator& gen) {
        std::uniform_real_distribution<float> dist{};
        return Scalar(dist(gen));
    }

    template<class RandomGenerator>
    inline Scalar<Float> Scalar<Float>::random_normal(RandomGenerator& gen) {
        std::normal_distribution<float> dist{};
        return Scalar(dist(gen));
    }

    template<class Distribution, class RandomGenerator>
    inline Scalar<Float> Scalar<Float>::random_any(Distribution& dist, RandomGenerator& gen) {
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

    template<class RandomGenerator>
    inline Scalar<Double> Scalar<Double>::random_uniform(RandomGenerator& gen) {
        std::uniform_real_distribution<double> dist{};
        return Scalar(dist(gen));
    }

    template<class RandomGenerator>
    inline Scalar<Double> Scalar<Double>::random_normal(RandomGenerator& gen) {
        std::normal_distribution<double> dist{};
        return Scalar(dist(gen));
    }

    template<class Distribution, class RandomGenerator>
    inline Scalar<Double> Scalar<Double>::random_any(Distribution& dist, RandomGenerator& gen) {
        return Scalar(dist(gen));
    }
}

#include "ScalarArithmetic.h"
#include "ElementaryFunction.h"
#include "Operation/Pow.h"

namespace Physica::Core {
    template<ScalarOption Option>
    inline Scalar<Option> operator^(const Scalar<Option>& s1, const Scalar<Option>& s2) {
        return pow(s1, s2);
    }

    template<>
    inline Scalar<MultiPrecision> operator^(const Scalar<MultiPrecision>& s1, const Scalar<MultiPrecision>& s2) {
        return s1.isInteger() ? powScalar(s1, s2) : exp(ln(s1) * s2);
    }
}
