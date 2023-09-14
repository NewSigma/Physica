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
    template<class ScalarType, DiffMode Mode>
    __host__ __device__ inline Differentiable<ScalarType, Mode> abs(const Differentiable<ScalarType, Mode>& s) {
        const ScalarType value = square(s.getValue());
        if constexpr (Mode == DiffMode::Forward)
            return {abs(s.getValue()), s.getValue().isPositive() ? s.getTangent() : -s.getTangent()};
        else {
            auto& tracer = DiffTracer<ScalarType>::getInstance();
            const size_t index = tracer.pushOperation(value, ExpressionType::Abs);
            tracer.pushOperand(s.getTraceIndex());
            return {value, index};
        }
    }

    template<class ScalarType, DiffMode Mode>
    __host__ __device__ inline Differentiable<ScalarType, Mode> square(const Differentiable<ScalarType, Mode>& s) {
        const ScalarType value = square(s.getValue());
        if constexpr (Mode == DiffMode::Forward)
            return {value, ScalarType(2) * s.getValue() * s.getTangent()};
        else {
            auto& tracer = DiffTracer<ScalarType>::getInstance();
            const size_t index = tracer.pushOperation(value, ExpressionType::Square);
            tracer.pushOperand(s.getTraceIndex());
            return {value, index};
        }
    }

    template<class ScalarType, DiffMode Mode>
    __host__ __device__ inline Differentiable<ScalarType, Mode> reciprocal(const Differentiable<ScalarType, Mode>& s) {
        const ScalarType rep = reciprocal(s.getValue());
        return {rep, -s.getTangent() * square(rep)};
    }

    template<class ScalarType, DiffMode Mode>
    __host__ __device__ Differentiable<ScalarType, Mode> sqrt(const Differentiable<ScalarType, Mode>& s) {
        const ScalarType value = sqrt(s.getValue());
        return {value, ScalarType(0.5) * s.getTangent() / value};
    }

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> cbrt(const Differentiable<ScalarType, Mode>& s) {
        constexpr double Factor = 1.0 / 3;
        const ScalarType value = cbrt(s.getValue());
        return {value, ScalarType(Factor) * value * s.getTangent() / s.getValue()};
    }

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> ln(const Differentiable<ScalarType, Mode>& s) {
        return {ln(s.getValue()), s.getTangent() / s.getValue()};
    }

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> exp(const Differentiable<ScalarType, Mode>& s) {
        const ScalarType value = exp(s.getValue());
        return {value, value * s.getTangent()};
    }

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> cos(const Differentiable<ScalarType, Mode>& s) {
        ScalarType sin_value, cos_value;
        sincos(s.getValue(), sin_value, cos_value);
        return {cos_value, -sin_value * s.getTangent()};
    }

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> sin(const Differentiable<ScalarType, Mode>& s) {
        ScalarType sin_value, cos_value;
        sincos(s.getValue(), sin_value, cos_value);
        return {sin_value, cos_value * s.getTangent()};
    }

    template<class ScalarType, DiffMode Mode>
    void sincos(Differentiable<ScalarType, Mode> s, Differentiable<ScalarType, Mode>& sin_result, Differentiable<ScalarType, Mode>& cos_result) {
        ScalarType sin_value, cos_value;
        sincos(s.getValue(), sin_value, cos_value);
        sin_result = Differentiable<ScalarType, Mode>(cos_value, -sin_value * s.getTangent());
        cos_result = Differentiable<ScalarType, Mode>(sin_value, cos_value * s.getTangent());
    }
}
