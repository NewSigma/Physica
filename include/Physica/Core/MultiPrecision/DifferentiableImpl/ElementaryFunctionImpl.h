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
    __host__ __device__ inline Differentiable<ScalarType, Mode, 1> abs(const Differentiable<ScalarType, Mode, 1>& s) {
        const ScalarType value = abs(s.getValue());
        if constexpr (Mode == DiffMode::Forward)
            return {value, s.getValue().isPositive() ? s.getGrad() : -s.getGrad()};
        else {
            auto& tracer = DiffTracer<ScalarType, 1>::getInstance();
            const auto result = tracer.pushOperation(value, ExpressionType::Abs);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode>
    inline Differentiable<ScalarType, Mode, 1> relu(const Differentiable<ScalarType, Mode, 1>& s) {
        const ScalarType value = relu(s.getValue());
        if constexpr (Mode == DiffMode::Forward)
            return {value, s.getValue().isPositive() ? s.getGrad() : ScalarType(0)};
        else {
            auto& tracer = DiffTracer<ScalarType, 1>::getInstance();
            const auto result = tracer.pushOperation(value, ExpressionType::Relu);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode>
    __host__ __device__ inline Differentiable<ScalarType, Mode, 1> square(const Differentiable<ScalarType, Mode, 1>& s) {
        const ScalarType value = square(s.getValue());
        if constexpr (Mode == DiffMode::Forward)
            return {value, ScalarType(2) * s.getValue() * s.getGrad()};
        else {
            auto& tracer = DiffTracer<ScalarType, 1>::getInstance();
            const auto result = tracer.pushOperation(value, ExpressionType::Square);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode>
    __host__ __device__ inline Differentiable<ScalarType, Mode, 1> reciprocal(const Differentiable<ScalarType, Mode, 1>& s) {
        const ScalarType rep = reciprocal(s.getValue());
        if constexpr (Mode == DiffMode::Forward)
            return {rep, -s.getGrad() * square(rep)};
        else {
            auto& tracer = DiffTracer<ScalarType, 1>::getInstance();
            const auto result = tracer.pushOperation(rep, ExpressionType::Reciprocal);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode>
    __host__ __device__ Differentiable<ScalarType, Mode, 1> sqrt(const Differentiable<ScalarType, Mode, 1>& s) {
        const ScalarType value = sqrt(s.getValue());
        if constexpr (Mode == DiffMode::Forward)
            return {value, ScalarType(0.5) * s.getGrad() / value};
        else {
            auto& tracer = DiffTracer<ScalarType, 1>::getInstance();
            const auto result = tracer.pushOperation(value, ExpressionType::Sqrt);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> cbrt(const Differentiable<ScalarType, Mode, 1>& s) {
        const ScalarType value = cbrt(s.getValue());
        if constexpr (Mode == DiffMode::Forward) {
            constexpr double Factor = 1.0 / 3;
            return {value, ScalarType(Factor) * value * s.getGrad() / s.getValue()};
        }
        else {
            auto& tracer = DiffTracer<ScalarType, 1>::getInstance();
            const auto result = tracer.pushOperation(value, ExpressionType::Cbrt);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> ln(const Differentiable<ScalarType, Mode, 1>& s) {
        const ScalarType value = ln(s.getValue());
        if constexpr (Mode == DiffMode::Forward)
            return {value, s.getGrad() / s.getValue()};
        else {
            auto& tracer = DiffTracer<ScalarType, 1>::getInstance();
            const auto result = tracer.pushOperation(value, ExpressionType::Ln);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> exp(const Differentiable<ScalarType, Mode, 1>& s) {
        const ScalarType value = exp(s.getValue());
        if constexpr (Mode == DiffMode::Forward)
            return {value, value * s.getGrad()};
        else {
            auto& tracer = DiffTracer<ScalarType, 1>::getInstance();
            const auto result = tracer.pushOperation(value, ExpressionType::Exp);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> cos(const Differentiable<ScalarType, Mode, 1>& s) {
        if constexpr (Mode == DiffMode::Forward) {
            ScalarType sin_value, cos_value;
            sincos(s.getValue(), sin_value, cos_value);
            return {cos_value, -sin_value * s.getGrad()};
        }
        else {
            const ScalarType value = cos(s.getValue());
            auto& tracer = DiffTracer<ScalarType, 1>::getInstance();
            const auto result = tracer.pushOperation(value, ExpressionType::Cos);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> sin(const Differentiable<ScalarType, Mode, 1>& s) {
        if constexpr (Mode == DiffMode::Forward) {
            ScalarType sin_value, cos_value;
            sincos(s.getValue(), sin_value, cos_value);
            return {sin_value, cos_value * s.getGrad()};
        }
        else {
            const ScalarType value = sin(s.getValue());
            auto& tracer = DiffTracer<ScalarType, 1>::getInstance();
            const auto result = tracer.pushOperation(value, ExpressionType::Sin);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode>
    void sincos(Differentiable<ScalarType, Mode, 1> s, Differentiable<ScalarType, Mode, 1>& sin_result, Differentiable<ScalarType, Mode, 1>& cos_result) {
        ScalarType sin_value, cos_value;
        sincos(s.getValue(), sin_value, cos_value);
        if constexpr (Mode == DiffMode::Forward) {
            sin_result = Differentiable<ScalarType, Mode, 1>(sin_value, cos_value * s.getGrad());
            cos_result = Differentiable<ScalarType, Mode, 1>(cos_value, -sin_value * s.getGrad());
        }
        else {
            auto& tracer = DiffTracer<ScalarType, 1>::getInstance();
            sin_result = tracer.pushOperation(sin_value, ExpressionType::Sin);
            tracer.pushOperand(s);
            cos_result = tracer.pushOperation(cos_value, ExpressionType::Cos);
            tracer.pushOperand(s);
        }
    }

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> arccos(const Differentiable<ScalarType, Mode, 1>& s) {
        if constexpr (Mode == DiffMode::Forward) {
            return {arccos(s.getValue()), -s.getGrad() / sqrt(ScalarType(1) - square(s.getValue()))};
        }
        else {
            const ScalarType value = arccos(s.getValue());
            auto& tracer = DiffTracer<ScalarType, 1>::getInstance();
            const auto result = tracer.pushOperation(value, ExpressionType::ArcCos);
            tracer.pushOperand(s);
            return result;
        }
    }
}
