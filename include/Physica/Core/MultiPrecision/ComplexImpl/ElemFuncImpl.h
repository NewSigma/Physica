/*
 * Copyright 2020-2021 Weibo He.
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
    ScalarType abs(const Complex<ScalarType>& c) {
        return c.norm();
    }

    template<class ScalarType>
    Complex<ScalarType> square(const Complex<ScalarType>& c) {
        const auto& real = c.getReal();
        const auto& imagine = c.getImag();
        return Complex<ScalarType>(square(real) - square(imagine), (real * imagine) << 1);
    }

    template<class ScalarType>
    inline Complex<ScalarType> reciprocal(const Complex<ScalarType>& c) {
        const auto& real = c.getReal();
        const auto& imagine = c.getImag();
        const auto divisor = reciprocal(square(real) + square(imagine));
        return Complex<ScalarType>(real * divisor, -imagine * divisor);
    }

    template<class ScalarType>
    Complex<ScalarType> sqrt(const Complex<ScalarType>& c) {
        using RealType = typename ScalarType::RealType;
        const RealType n = sqrt(c.norm());
        const RealType a = c.phase() * RealType(0.5);
        return Complex<ScalarType>(n * cos(a), n * sin(a));
    }

    template<class ScalarType>
    Complex<ScalarType> ln(const Complex<ScalarType>& c) {
        const auto n = c.norm();
        const auto a = c.phase();
        return Complex<ScalarType>(ln(n), a);
    }

    template<class ScalarType>
    Complex<ScalarType> exp(const Complex<ScalarType>& c) {
        const auto& exp_real = exp(c.getReal());
        const auto& imagine = c.getImag();
        return Complex<ScalarType>(exp_real * cos(imagine), exp_real * sin(imagine));
    }

    template<class ScalarType>
    Complex<ScalarType> cos(const Complex<ScalarType>& c) {
        const auto& real = c.getReal();
        const auto& imagine = c.getImag();
        return Complex<ScalarType>(cos(real) * cosh(imagine), - sin(real) * sinh(imagine));
    }

    template<class ScalarType>
    Complex<ScalarType> sin(const Complex<ScalarType>& c) {
        const auto& real = c.getReal();
        const auto& imagine = c.getImag();
        return Complex<ScalarType>(sin(real) * cosh(imagine), cos(real) * sinh(imagine));
    }

    template<class ScalarType>
    Complex<ScalarType> tan(const Complex<ScalarType>& c) {
        return sin(c) / cos(c);
    }

    template<class ScalarType>
    Complex<ScalarType> sec(const Complex<ScalarType>& c) {
        return reciprocal(cos(c));
    }

    template<class ScalarType>
    Complex<ScalarType> csc(const Complex<ScalarType>& c) {
        return reciprocal(sin(c));
    }

    template<class ScalarType>
    Complex<ScalarType> cot(const Complex<ScalarType>& c) {
        return reciprocal(tan(c));
    }

    template<class ScalarType>
    Complex<ScalarType> cosh(const Complex<ScalarType>& c) {
        return (exp(c) + exp(-c)) >> 1;
    }

    template<class ScalarType>
    Complex<ScalarType> sinh(const Complex<ScalarType>& c) {
        return (exp(c) - exp(-c)) >> 1;
    }

    template<class ScalarType>
    Complex<ScalarType> tanh(const Complex<ScalarType>& c) {
        return sinh(c) / cosh(c);
    }

    template<class ScalarType>
    Complex<ScalarType> sech(const Complex<ScalarType>& c) {
        return reciprocal(cosh(c));
    }

    template<class ScalarType>
    Complex<ScalarType> csch(const Complex<ScalarType>& c) {
        return reciprocal(sinh(c));
    }

    template<class ScalarType>
    Complex<ScalarType> coth(const Complex<ScalarType>& c) {
        return reciprocal(tanh(c));
    }
}
