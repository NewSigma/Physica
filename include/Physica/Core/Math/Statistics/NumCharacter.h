/*
 * Copyright 2022-2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"

namespace Physica::Core {
    template<Vector T>
    T::ScalarType mean(const T& x) {
        using ScalarType = T::ScalarType;
        return x.sum() / ScalarType(x.getLength());
    }

    template<Scalar T>
    __host__ __device__ inline void toNextMean(T& mean, size_t lastNumSample, T sample) {
        const T factor1 = T(lastNumSample);
        const T factor2 = reciprocal(T(lastNumSample + 1));
        mean = (factor1 * mean + sample) * factor2;
    }

    template<Vector T1, Vector T2>
    inline void toNextMean(T1& mean, size_t lastNumSample, const T2& sample) {
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename T1::ScalarType, typename T2::ScalarType>::Type;
        const ScalarType factor1 = ScalarType(lastNumSample);
        const ScalarType factor2 = reciprocal(ScalarType(lastNumSample + 1));
        mean = (factor1 * mean + sample) * factor2;
    }

    template<Matrix T1, Matrix T2>
    inline void toNextMean(T1& mean, size_t lastNumSample, const T2& sample) {
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename T1::ScalarType, typename T2::ScalarType>::Type;
        const ScalarType factor1 = ScalarType(lastNumSample);
        const ScalarType factor2 = reciprocal(ScalarType(lastNumSample + 1));
        mean = (factor1 * mean + sample) * factor2;
    }

    template<Vector T>
    T::ScalarType mean_stable(const T& x) {
        using ScalarType = T::ScalarType;
        ScalarType result = 0;
        for (size_t i = 0; i < x.getLength(); ++i)
            toNextMean(result, i, x.calc(i));
        return result;
    }

    template<Vector T>
    T::ScalarType variance(const T& x) {
        using ScalarType = T::ScalarType;
        const size_t length = x.getLength();
        const ScalarType x_mean = mean(x);
        return square(x - x_mean).sum() / ScalarType(length - 1);
    }

    template<Vector T>
    T::ScalarType variance(const T& x, typename T::ScalarType prior_mean) {
        using ScalarType = T::ScalarType;
        const size_t length = x.getLength();
        return square(x - prior_mean).sum() / ScalarType(length);
    }

    template<Scalar T>
    inline void toNextVariance(T& var, T& mean, size_t lastNumSample, T sample) {
        const T factor1 = T(lastNumSample);
        const T factor2 = reciprocal(T(lastNumSample + 1));
        var = (var + square(mean - sample) * factor2) * (factor1 * factor2);
        toNextMean(mean, lastNumSample, sample);
    }

    template<Vector T1, Vector T2>
    inline void toNextVariance(T1& var, T1& mean, size_t lastNumSample, const T2& sample) {
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename T1::ScalarType, typename T2::ScalarType>::Type;
        const ScalarType factor1 = ScalarType(lastNumSample);
        const ScalarType factor2 = reciprocal(ScalarType(lastNumSample + 1));
        var = (var + square(mean - sample) * factor2) * (factor1 * factor2);
        toNextMean(mean, lastNumSample, sample);
    }

    template<Matrix T1, Matrix T2>
    inline void toNextVariance(T1& var, T1& mean, size_t lastNumSample, const T2& sample) {
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename T1::ScalarType, typename T2::ScalarType>::Type;
        const ScalarType factor1 = ScalarType(lastNumSample);
        const ScalarType factor2 = reciprocal(ScalarType(lastNumSample + 1));
        var = (var + square_elem(mean - sample) * factor2) * (factor1 * factor2);
        toNextMean(mean, lastNumSample, sample);
    }
    /**
     * Stable if large dataset is used, prior version is not provided because they behave similarly at large dataset.
     */
    template<Vector T>
    T::ScalarType variance_stable(const T& x) {
        using ScalarType = T::ScalarType;
        ScalarType result = 0;
        ScalarType mean = 0;
        for (size_t i = 0; i < x.getLength(); ++i)
            toNextVariance(result, mean, i, x[i]);
        return result;
    }

    template<Vector T>
    inline T::ScalarType deviation(const T& x) {
        return sqrt(variance(x));
    }

    template<Vector T>
    inline T::ScalarType deviation_stable(const T& x) {
        return sqrt(variance_stable(x));
    }

    template<Vector T>
    DenseVector<typename T::ScalarType, T::SizeAtCompile>
    normalize(const T& x) {
        using ScalarType = T::ScalarType;
        const ScalarType x_mean = mean(x);
        const ScalarType factor = reciprocal(deviation(x));
        return (x - x_mean) * factor;
    }

    template<Vector T>
    T::ScalarType covariance(const T& x, const T& y) {
        assert(x.getLength() == y.getLength());
        using ScalarType = T::ScalarType;
        const ScalarType x_mean = mean(x);
        const ScalarType y_mean = mean(y);
        DenseVector<typename T::ScalarType, T::SizeAtCompile> temp = hadamard((x - x_mean), (y - y_mean));
        return temp.sum() / ScalarType(temp.getLength() - 1);
    }

    template<Vector T>
    T::ScalarType skew(const T& x) {
        using ScalarType = T::ScalarType;
        T temp = normalize(x);
        temp = hadamard(square(temp), temp);
        const size_t length = x.getLength();
        const ScalarType factor = ScalarType(length * length) / ScalarType((length - 1) * (length - 2));
        return mean(temp) * factor;
    }

    template<Vector T>
    T::ScalarType excess_kurt(const T& x) {
        using ScalarType = T::ScalarType;
        T temp = normalize(x);
        temp = square(temp);
        const ScalarType mean2 = mean(temp);
        temp = square(temp);
        const ScalarType mean1 = mean(temp);

        const size_t length = x.getLength();
        const ScalarType factor2 = ScalarType(length * length * 3) / ScalarType((length - 2) * (length - 3));
        const ScalarType factor1 = ScalarType(length * length * (length + 1)) / ScalarType((length - 1) * (length - 2) * (length - 3));
        return factor1 * mean1 - factor2 * mean2;
    }

    template<Vector T>
    inline T::ScalarType kurt(const T& x) {
        return excess_kurt(x) + 3.0;
    }
}
