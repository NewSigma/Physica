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
    template<class ScalarType, DiffMode Mode, class OtherScalar>
    [[nodiscard]] inline typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type
    operator+(const Differentiable<ScalarType, Mode>& s1, const ScalarBase<OtherScalar>& s2) {
        using ResultType = typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type;
        if constexpr (OtherScalar::isDifferentiable)
            return ResultType(s1.getValue() + s2.getValue(), s1.getTangent() + s2.getTangent());
        else
            return ResultType(s1.getValue() + s2.getValue(), s1.getTangent());
    }

    template<class ScalarType, DiffMode Mode, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type>::type
    operator+(const ScalarBase<OtherScalar>& s1, const Differentiable<ScalarType, Mode>& s2) {
        return s2 + s1.getDerived();
    }

    template<class ScalarType, DiffMode Mode, class OtherScalar>
    [[nodiscard]] inline typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type
    operator-(const Differentiable<ScalarType, Mode>& s1, const ScalarBase<OtherScalar>& s2) {
        using ResultType = typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type;
        if constexpr (OtherScalar::isDifferentiable)
            return ResultType(s1.getValue() - s2.getValue(), s1.getTangent() - s2.getTangent());
        else
            return ResultType(s1.getValue() - s2.getValue(), s1.getTangent());
    }

    template<class ScalarType, DiffMode Mode, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type>::type
    operator-(const ScalarBase<OtherScalar>& s1, const Differentiable<ScalarType, Mode>& s2) {
        return -(s2 - s1.getDerived());
    }

    template<class ScalarType, DiffMode Mode, class OtherScalar>
    [[nodiscard]] inline typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type
    operator*(const Differentiable<ScalarType, Mode>& s1, const ScalarBase<OtherScalar>& s2) {
        using ResultType = typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type;
        if constexpr (OtherScalar::isDifferentiable)
            return ResultType(s1.getValue() * s2.getValue(), s1.getTangent() * s2.getValue() + s1.getValue() * s2.getTangent());
        else
            return ResultType(s1.getValue() * s2.getValue(), s1.getTangent() * s2.getValue());
    }

    template<class ScalarType, DiffMode Mode, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type>::type
    operator*(const ScalarBase<OtherScalar>& s1, const Differentiable<ScalarType, Mode>& s2) {
        return s2 * s1.getDerived();
    }

    template<class ScalarType, DiffMode Mode, class OtherScalar>
    [[nodiscard]] inline typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type
    operator/(const Differentiable<ScalarType, Mode>& s1, const ScalarBase<OtherScalar>& s2) {
        using ResultType = typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type;
        const auto rep = reciprocal(s2.getValue());
        if constexpr (OtherScalar::isDifferentiable)
            return ResultType(s1.getValue() * rep, (s1.getTangent() * s2.getValue() - s1.getValue() * s2.getTangent()) * square(rep));
        else
            return ResultType(s1.getValue() * rep, s1.getTangent() * rep);
    }

    template<class ScalarType, DiffMode Mode, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type>::type
    operator/(const ScalarBase<OtherScalar>& s1, const Differentiable<ScalarType, Mode>& s2) {
        using ResultType = typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type;
        const auto rep = reciprocal(s2.getValue());
        return ResultType(s1.getValue() * rep, -s1.getValue() * s2.getTangent() * square(rep));
    }
}
