/*
 * Copyright 2020-2021 WeiBo He.
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
    template<class ScalarType>
    ScalarType abs(const ComplexScalar<ScalarType>& c) {
        return c.norm();
    }

    template<class ScalarType>
    ComplexScalar<ScalarType> square(const ComplexScalar<ScalarType>& c) {
        const auto& real = c.getReal();
        const auto& imagine = c.getImag();
        return ComplexScalar<ScalarType>(square(real) - square(imagine), (real * imagine) << 1);
    }

    template<class ScalarType>
    inline ComplexScalar<ScalarType> reciprocal(const ComplexScalar<ScalarType>& c) {
        const auto& real = c.getReal();
        const auto& imagine = c.getImag();
        const auto divisor = reciprocal(square(real) + square(imagine));
        return ComplexScalar<ScalarType>(real * divisor, -imagine * divisor);
    }

    template<class ScalarType>
    ComplexScalar<ScalarType> sqrt(const ComplexScalar<ScalarType>& c) {
        using RealType = typename ScalarType::RealType;
        const RealType n = sqrt(c.norm());
        const RealType a = c.arg() * RealType(0.5);
        return ComplexScalar<ScalarType>(n * cos(a), n * sin(a));
    }

    template<class ScalarType>
    ComplexScalar<ScalarType> ln(const ComplexScalar<ScalarType>& c) {
        const auto n = norm(c);
        const auto a = arg(c);
        return ComplexScalar<ScalarType>(ln(n), a);
    }

    template<class ScalarType>
    ComplexScalar<ScalarType> exp(const ComplexScalar<ScalarType>& c) {
        const auto& exp_real = exp(c.getReal());
        const auto& imagine = c.getImag();
        return ComplexScalar<ScalarType>(exp_real * cos(imagine), exp_real * sin(imagine));
    }

    template<class ScalarType>
    ComplexScalar<ScalarType> cos(const ComplexScalar<ScalarType>& c) {
        const auto& real = c.getReal();
        const auto& imagine = c.getImag();
        return ComplexScalar<ScalarType>(cos(real) * cosh(imagine), - sin(real) * sinh(imagine));
    }

    template<class ScalarType>
    ComplexScalar<ScalarType> sin(const ComplexScalar<ScalarType>& c) {
        const auto& real = c.getReal();
        const auto& imagine = c.getImag();
        return ComplexScalar<ScalarType>(sin(real) * cosh(imagine), cos(real) * sinh(imagine));
    }

    template<class ScalarType>
    ComplexScalar<ScalarType> tan(const ComplexScalar<ScalarType>& c) {
        return sin(c) / cos(c);
    }

    template<class ScalarType>
    ComplexScalar<ScalarType> sec(const ComplexScalar<ScalarType>& c) {
        return reciprocal(cos(c));
    }

    template<class ScalarType>
    ComplexScalar<ScalarType> csc(const ComplexScalar<ScalarType>& c) {
        return reciprocal(sin(c));
    }

    template<class ScalarType>
    ComplexScalar<ScalarType> cot(const ComplexScalar<ScalarType>& c) {
        return reciprocal(tan(c));
    }

    template<class ScalarType>
    ComplexScalar<ScalarType> cosh(const ComplexScalar<ScalarType>& c) {
        return (exp(c) + exp(-c)) >> 1;
    }

    template<class ScalarType>
    ComplexScalar<ScalarType> sinh(const ComplexScalar<ScalarType>& c) {
        return (exp(c) - exp(-c)) >> 1;
    }

    template<class ScalarType>
    ComplexScalar<ScalarType> tanh(const ComplexScalar<ScalarType>& c) {
        return sinh(c) / cosh(c);
    }

    template<class ScalarType>
    ComplexScalar<ScalarType> sech(const ComplexScalar<ScalarType>& c) {
        return reciprocal(cosh(c));
    }

    template<class ScalarType>
    ComplexScalar<ScalarType> csch(const ComplexScalar<ScalarType>& c) {
        return reciprocal(sinh(c));
    }

    template<class ScalarType>
    ComplexScalar<ScalarType> coth(const ComplexScalar<ScalarType>& c) {
        return reciprocal(tanh(c));
    }
}
