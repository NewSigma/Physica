/*
 * Copyright 2023 WeiBo He.
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

namespace Physica::Core {
    template<class ScalarType, class OtherScalar>
    [[nodiscard]] inline Differentiable<typename Internal::BinaryScalarOpReturnType<ScalarType, OtherScalar>::Type>
    operator+(const Differentiable<ScalarType>& s1, const OtherScalar& s2) {
        if constexpr (OtherScalar::isDifferentiable)
            return {s1.getValue() + s2.getValue(), s1.getTangent() + s2.getTangent()};
        else
            return {s1.getValue() + s2.getValue(), s1.getTangent()};
    }

    template<class ScalarType, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, Differentiable<typename Internal::BinaryScalarOpReturnType<ScalarType, OtherScalar>::Type>>::type
    operator+(const OtherScalar& s1, const Differentiable<ScalarType>& s2) {
        return s2 + s1;
    }

    template<class ScalarType, class OtherScalar>
    [[nodiscard]] inline Differentiable<typename Internal::BinaryScalarOpReturnType<ScalarType, OtherScalar>::Type>
    operator-(const Differentiable<ScalarType>& s1, const OtherScalar& s2) {
        if constexpr (OtherScalar::isDifferentiable)
            return {s1.getValue() - s2.getValue(), s1.getTangent() - s2.getTangent()};
        else
            return {s1.getValue() - s2.getValue(), s1.getTangent()};
    }

    template<class ScalarType, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, Differentiable<typename Internal::BinaryScalarOpReturnType<ScalarType, OtherScalar>::Type>>::type
    operator-(const OtherScalar& s1, const Differentiable<ScalarType>& s2) {
        return -(s2 - s1);
    }

    template<class ScalarType, class OtherScalar>
    [[nodiscard]] inline Differentiable<typename Internal::BinaryScalarOpReturnType<ScalarType, OtherScalar>::Type>
    operator*(const Differentiable<ScalarType>& s1, const OtherScalar& s2) {
        if constexpr (OtherScalar::isDifferentiable)
            return {s1.getValue() * s2.getValue(), s1.getTangent() * s2.getValue()};
        else
            return {s1.getValue() * s2.getValue(), s1.getTangent() * s2.getValue() + s1.getValue() * s2.getTangent()};
    }

    template<class ScalarType, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, Differentiable<typename Internal::BinaryScalarOpReturnType<ScalarType, OtherScalar>::Type>>::type
    operator*(const OtherScalar& s1, const Differentiable<ScalarType>& s2) {
        return s2 * s1;
    }

    template<class ScalarType, class OtherScalar>
    [[nodiscard]] inline Differentiable<typename Internal::BinaryScalarOpReturnType<ScalarType, OtherScalar>::Type>
    operator/(const Differentiable<ScalarType>& s1, const OtherScalar& s2) {
        using ResultType = typename Internal::BinaryScalarOpReturnType<ScalarType, OtherScalar>::Type;
        const ResultType rep = reciprocal(s2.getValue());
        if constexpr (OtherScalar::isDifferentiable)
            return {s1.getValue() * rep, s1.getTangent() * rep};
        else
            return {s1.getValue() * rep, (s1.getTangent() * s2.getValue() - s1.getValue() * s2.getTangent()) * square(rep)};
    }

    template<class ScalarType, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, Differentiable<typename Internal::BinaryScalarOpReturnType<ScalarType, OtherScalar>::Type>>::type
    operator/(const OtherScalar& s1, const Differentiable<ScalarType>& s2) {
        using ResultType = typename Internal::BinaryScalarOpReturnType<ScalarType, OtherScalar>::Type;
        const ResultType rep = reciprocal(s2.getValue());
        return {s1.getValue() * rep, -s1.getValue() * s2.getTangent() * square(rep)};
    }

    template<class ScalarType>
    __host__ __device__ inline Differentiable<ScalarType> abs(const Differentiable<ScalarType>& s) {
        return {abs(s.getValue()), s.getValue().isPositive() ? s.getTangent() : -s.getTangent()};
    }

    template<class ScalarType>
    __host__ __device__ inline Differentiable<ScalarType> square(const Differentiable<ScalarType>& s) {
        return {square(s.getValue()), ScalarType(2) * s.getValue() * s.getTangent()};
    }

    template<class ScalarType>
    __host__ __device__ inline Differentiable<ScalarType> reciprocal(const Differentiable<ScalarType>& s) {
        const ScalarType rep = reciprocal(s2.getValue());
        return {rep, -s.getTangent() * square(rep)};
    }

    template<class ScalarType>
    __host__ __device__ Differentiable<ScalarType> sqrt(const Differentiable<ScalarType>& s) {
        const ScalarType value = sqrt(s.getValue());
        return {value, ScalarType(0.5) * s.getTangent() / value};
    }

    template<class ScalarType>
    Differentiable<ScalarType> cbrt(const Differentiable<ScalarType>& s) {
        constexpr double Factor = 1.0 / 3;
        const ScalarType value = cbrt(s.getValue());
        return {value, ScalarType(Factor) * value * s.getTangent() / s.getValue()};
    }

    template<class ScalarType>
    Differentiable<ScalarType> ln(const Differentiable<ScalarType>& s) {
        return {ln(s.getValue()), s.getTangent() / s.getValue()};
    }

    template<class ScalarType>
    Differentiable<ScalarType> exp(const Differentiable<ScalarType>& s) {
        const ScalarType value = exp(s.getValue());
        return {value, value * s.getTangent()};
    }

    template<class ScalarType>
    Differentiable<ScalarType> cos(const Differentiable<ScalarType>& s) {
        ScalarType sin_value, cos_value;
        sincos(s.getValue(), sin_value, cos_value);
        return {cos_value, -sin_value * s.getTangent()};
    }

    template<class ScalarType>
    Differentiable<ScalarType> sin(const Differentiable<ScalarType>& s) {
        ScalarType sin_value, cos_value;
        sincos(s.getValue(), sin_value, cos_value);
        return {sin_value, cos_value * s.getTangent()};
    }

    template<class ScalarType>
    void sincos(Differentiable<ScalarType> s, Differentiable<ScalarType>& sin_result, Differentiable<ScalarType>& cos_result) {
        ScalarType sin_value, cos_value;
        sincos(s.getValue(), sin_value, cos_value);
        sin_result = Differentiable<ScalarType>(cos_value, -sin_value * s.getTangent());
        cos_result = Differentiable<ScalarType>(sin_value, cos_value * s.getTangent());
    }
}
