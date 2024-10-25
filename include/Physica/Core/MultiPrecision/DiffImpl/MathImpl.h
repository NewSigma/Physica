/*
 * Copyright 2023-2024 Weibo He.
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
    template<class ScalarType, DiffMode Mode, unsigned int Order>
    __host__ __device__ inline Diff<ScalarType, Mode, Order> abs(const Diff<ScalarType, Mode, Order>& s) {
        const ScalarType value = abs(s.getValue());
        if constexpr (Mode == DiffMode::Forward)
            return {value, s.getValue().isPositive() ? s.getGrad() : -s.getGrad()};
        else {
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            const auto result = tracer.pushOperation(value, ExprType::Abs);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    inline Diff<ScalarType, Mode, Order> relu(const Diff<ScalarType, Mode, Order>& s) {
        const ScalarType value = relu(s.getValue());
        if constexpr (Mode == DiffMode::Forward)
            return {value, s.getValue().isPositive() ? s.getGrad() : ScalarType(0)};
        else {
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            const auto result = tracer.pushOperation(value, ExprType::Relu);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    __host__ __device__ inline Diff<ScalarType, Mode, Order> square(const Diff<ScalarType, Mode, Order>& s) {
        const ScalarType value = square(s.getValue());
        if constexpr (Mode == DiffMode::Forward)
            return {value, ScalarType(2) * s.getValue() * s.getGrad()};
        else {
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            const auto result = tracer.pushOperation(value, ExprType::Square);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    __host__ __device__ inline Diff<ScalarType, Mode, Order> reciprocal(const Diff<ScalarType, Mode, Order>& s) {
        const ScalarType rep = reciprocal(s.getValue());
        if constexpr (Mode == DiffMode::Forward)
            return {rep, -s.getGrad() * square(rep)};
        else {
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            const auto result = tracer.pushOperation(rep, ExprType::Reciprocal);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    __host__ __device__ Diff<ScalarType, Mode, Order> sqrt(const Diff<ScalarType, Mode, Order>& s) {
        const ScalarType value = sqrt(s.getValue());
        if constexpr (Mode == DiffMode::Forward)
            return {value, ScalarType(0.5) * s.getGrad() / value};
        else {
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            const auto result = tracer.pushOperation(value, ExprType::Sqrt);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Diff<ScalarType, Mode, Order> cbrt(const Diff<ScalarType, Mode, Order>& s) {
        const ScalarType value = cbrt(s.getValue());
        if constexpr (Mode == DiffMode::Forward) {
            constexpr double Factor = 1.0 / 3;
            return {value, ScalarType(Factor) * value * s.getGrad() / s.getValue()};
        }
        else {
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            const auto result = tracer.pushOperation(value, ExprType::Cbrt);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Diff<ScalarType, Mode, Order> ln(const Diff<ScalarType, Mode, Order>& s) {
        const ScalarType value = ln(s.getValue());
        if constexpr (Mode == DiffMode::Forward)
            return {value, s.getGrad() / s.getValue()};
        else {
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            const auto result = tracer.pushOperation(value, ExprType::Ln);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Diff<ScalarType, Mode, Order> exp(const Diff<ScalarType, Mode, Order>& s) {
        const ScalarType value = exp(s.getValue());
        if constexpr (Mode == DiffMode::Forward)
            return {value, value * s.getGrad()};
        else {
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            const auto result = tracer.pushOperation(value, ExprType::Exp);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Diff<ScalarType, Mode, Order> cos(const Diff<ScalarType, Mode, Order>& s) {
        if constexpr (Mode == DiffMode::Forward) {
            ScalarType sin_value, cos_value;
            sincos(s.getValue(), sin_value, cos_value);
            return {cos_value, -sin_value * s.getGrad()};
        }
        else {
            const ScalarType value = cos(s.getValue());
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            const auto result = tracer.pushOperation(value, ExprType::Cos);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Diff<ScalarType, Mode, Order> sin(const Diff<ScalarType, Mode, Order>& s) {
        if constexpr (Mode == DiffMode::Forward) {
            ScalarType sin_value, cos_value;
            sincos(s.getValue(), sin_value, cos_value);
            return {sin_value, cos_value * s.getGrad()};
        }
        else {
            const ScalarType value = sin(s.getValue());
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            const auto result = tracer.pushOperation(value, ExprType::Sin);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    void sincos(Diff<ScalarType, Mode, Order> s, Diff<ScalarType, Mode, Order>& sin_result, Diff<ScalarType, Mode, Order>& cos_result) {
        ScalarType sin_value, cos_value;
        sincos(s.getValue(), sin_value, cos_value);
        if constexpr (Mode == DiffMode::Forward) {
            sin_result = Diff<ScalarType, Mode, Order>(sin_value, cos_value * s.getGrad());
            cos_result = Diff<ScalarType, Mode, Order>(cos_value, -sin_value * s.getGrad());
        }
        else {
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            sin_result = tracer.pushOperation(sin_value, ExprType::Sin);
            tracer.pushOperand(s);
            cos_result = tracer.pushOperation(cos_value, ExprType::Cos);
            tracer.pushOperand(s);
        }
    }

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Diff<ScalarType, Mode, Order> arccos(const Diff<ScalarType, Mode, Order>& s) {
        if constexpr (Mode == DiffMode::Forward) {
            return {arccos(s.getValue()), -s.getGrad() / sqrt(ScalarType(1) - square(s.getValue()))};
        }
        else {
            const ScalarType value = arccos(s.getValue());
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            const auto result = tracer.pushOperation(value, ExprType::ArcCos);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Diff<ScalarType, Mode, Order> tanh(const Diff<ScalarType, Mode, Order>& s) {
        const ScalarType value = tanh(s.getValue());
        if constexpr (Mode == DiffMode::Forward) {
            return {value, ScalarType(1) - square(value)};
        }
        else {
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            const auto result = tracer.pushOperation(value, ExprType::Tanh);
            tracer.pushOperand(s);
            return result;
        }
    }
}
