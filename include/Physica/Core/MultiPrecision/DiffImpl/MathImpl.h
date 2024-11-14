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
    template<class ScalarType, DiffMode Mode, int Order>
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

    template<class ScalarType, DiffMode Mode, int Order>
    __host__ __device__ inline auto abs(const ScalarRef<Diff<ScalarType, Mode, Order>>& x) {
        return abs(Diff<ScalarType, Mode, Order>(x));
    }

    template<class ScalarType, DiffMode Mode, int Order>
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

    template<class ScalarType, DiffMode Mode, int Order>
    __host__ __device__ inline Diff<ScalarType, Mode, Order> square(const Diff<ScalarType, Mode, Order>& s) {
        using ResultType = Diff<ScalarType, Mode, Order>;
        const ScalarType value = square(s.getValue());
        if constexpr (Mode == DiffMode::Forward) {
            using GradType = typename ResultType::GradType;
            return ResultType(value, GradType(ScalarType(2) * s * s.getGrad()));
        }
        else {
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            const auto result = tracer.pushOperation(value, ExprType::Square);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, int Order>
    __host__ __device__ inline auto square(const ScalarRef<Diff<ScalarType, Mode, Order>>& x) {
        return square(Diff<ScalarType, Mode, Order>(x));
    }

    template<class ScalarType, DiffMode Mode, int Order>
    __host__ __device__ inline Diff<ScalarType, Mode, Order> reciprocal(const Diff<ScalarType, Mode, Order>& s) {
        using ResultType = Diff<ScalarType, Mode, Order>;
        if constexpr (Mode == DiffMode::Forward) {
            using GradType = typename ResultType::GradType;
            const GradType v = reciprocal(GradType(s));
            return ResultType(v.getValue(), -s.getGrad() * square(v));
        }
        else {
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            const auto result = tracer.pushOperation(reciprocal(s.getValue()), ExprType::Reciprocal);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, int Order>
    __host__ __device__ Diff<ScalarType, Mode, Order> sqrt(const Diff<ScalarType, Mode, Order>& s) {
        using ResultType = Diff<ScalarType, Mode, Order>;
        if constexpr (Mode == DiffMode::Forward) {
            using GradType = typename ResultType::GradType;
            const GradType v = sqrt(GradType(s));
            return ResultType(v.getValue(), ScalarType(0.5) * s.getGrad() / v);
        }
        else {
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            const auto result = tracer.pushOperation(sqrt(s.getValue()), ExprType::Sqrt);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> cbrt(const Diff<ScalarType, Mode, Order>& s) {
        using ResultType = Diff<ScalarType, Mode, Order>;
        if constexpr (Mode == DiffMode::Forward) {
            using GradType = typename ResultType::GradType;
            const GradType v = cbrt(GradType(s));
            return ResultType(v.getValue(), ScalarType(1.0 / 3) * v * s.getGrad() / GradType(s));
        }
        else {
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            const auto result = tracer.pushOperation(cbrt(s.getValue()), ExprType::Cbrt);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> ln(const Diff<ScalarType, Mode, Order>& s) {
        using ResultType = Diff<ScalarType, Mode, Order>;
        if constexpr (Mode == DiffMode::Forward) {
            using GradType = typename ResultType::GradType;
            return ResultType(ln(s.getValue()), s.getGrad() / GradType(s));
        }
        else {
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            const auto result = tracer.pushOperation(ln(s.getValue()), ExprType::Ln);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> ln1p(const Diff<ScalarType, Mode, Order>& x) {
        static_assert(Mode == DiffMode::Forward, "[Error]: Not implemented");
        using ResultType = Diff<ScalarType, Mode, Order>;
        if constexpr (Mode == DiffMode::Forward) {
            using GradType = typename ResultType::GradType;
            return ResultType(ln1p(x.getValue()), x.getGrad() / (ScalarType(1) + GradType(x)));
        }
        else
            noImpl();
    }

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> exp(const Diff<ScalarType, Mode, Order>& s) {
        using ResultType = Diff<ScalarType, Mode, Order>;
        if constexpr (Mode == DiffMode::Forward) {
            using GradType = typename ResultType::GradType;
            const GradType v = exp(GradType(s));
            return {v.getValue(), v * s.getGrad()};
        }
        else {
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            const auto result = tracer.pushOperation(exp(s.getValue()), ExprType::Exp);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> cos(const Diff<ScalarType, Mode, Order>& s) {
        using ResultType = Diff<ScalarType, Mode, Order>;
        if constexpr (Mode == DiffMode::Forward) {
            using GradType = typename ResultType::GradType;
            GradType sin_value, cos_value;
            sincos(GradType(s), sin_value, cos_value);
            return ResultType(cos_value.getValue(), -sin_value * s.getGrad());
        }
        else {
            const ScalarType value = cos(s.getValue());
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            const auto result = tracer.pushOperation(value, ExprType::Cos);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> sin(const Diff<ScalarType, Mode, Order>& s) {
        using ResultType = Diff<ScalarType, Mode, Order>;
        if constexpr (Mode == DiffMode::Forward) {
            using GradType = typename ResultType::GradType;
            GradType sin_value, cos_value;
            sincos(GradType(s), sin_value, cos_value);
            return ResultType(sin_value.getValue(), cos_value * s.getGrad());
        }
        else {
            const ScalarType value = sin(s.getValue());
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            const auto result = tracer.pushOperation(value, ExprType::Sin);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, int Order>
    void sincos(Diff<ScalarType, Mode, Order> s, Diff<ScalarType, Mode, Order>& sin_result, Diff<ScalarType, Mode, Order>& cos_result) {
        using ResultType = Diff<ScalarType, Mode, Order>;
        if constexpr (Mode == DiffMode::Forward) {
            using GradType = typename ResultType::GradType;
            GradType sin_value, cos_value;
            sincos(GradType(s), sin_value, cos_value);
            sin_result = ResultType(sin_value, cos_value * s.getGrad());
            cos_result = ResultType(cos_value, -sin_value * s.getGrad());
        }
        else {
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            ScalarType sin_value, cos_value;
            sincos(s.getValue(), sin_value, cos_value);
            sin_result = tracer.pushOperation(sin_value, ExprType::Sin);
            tracer.pushOperand(s);
            cos_result = tracer.pushOperation(cos_value, ExprType::Cos);
            tracer.pushOperand(s);
        }
    }

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> tan(const Diff<ScalarType, Mode, Order>& s) {
        using ResultType = Diff<ScalarType, Mode, Order>;
        if constexpr (Mode == DiffMode::Forward) {
            using GradType = typename ResultType::GradType;
            return ResultType(tan(s.getValue()), s.getGrad() * square(sec(GradType(s))));
        }
        else
            noImpl();
    }

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> sec(const Diff<ScalarType, Mode, Order>& s) {
        using ResultType = Diff<ScalarType, Mode, Order>;
        if constexpr (Mode == DiffMode::Forward) {
            using GradType = typename ResultType::GradType;
            const auto x1 = GradType(s);
            const auto v = sec(x1);
            return ResultType(v.getValue(), s.getGrad() * v * tan(x1));
        }
        else
            noImpl();
    }

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> csc(const Diff<ScalarType, Mode, Order>& s) {
        using ResultType = Diff<ScalarType, Mode, Order>;
        if constexpr (Mode == DiffMode::Forward) {
            using GradType = typename ResultType::GradType;
            const auto x1 = GradType(s);
            const auto v = csc(x1);
            return ResultType(v.getValue(), -s.getGrad() * v * cot(x1));
        }
        else
            noImpl();
    }

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> cot(const Diff<ScalarType, Mode, Order>& s) {
        using ResultType = Diff<ScalarType, Mode, Order>;
        if constexpr (Mode == DiffMode::Forward) {
            using GradType = typename ResultType::GradType;
            return ResultType(cot(s.getValue()), -s.getGrad() * square(csc(GradType(s))));
        }
        else
            noImpl();
    }

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> arccos(const Diff<ScalarType, Mode, Order>& s) {
        using ResultType = Diff<ScalarType, Mode, Order>;
        if constexpr (Mode == DiffMode::Forward) {
            using GradType = typename ResultType::GradType;
            return ResultType(arccos(s.getValue()), -s.getGrad() / sqrt(ScalarType(1) - square(GradType(s))));
        }
        else {
            const ScalarType value = arccos(s.getValue());
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            const auto result = tracer.pushOperation(value, ExprType::ArcCos);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> cosh(const Diff<ScalarType, Mode, Order>& s) {
        using ResultType = Diff<ScalarType, Mode, Order>;
        if constexpr (Mode == DiffMode::Forward) {
            using GradType = typename ResultType::GradType;
            return ResultType(cosh(s.getValue()), s.getGrad() * sinh(GradType(s)));
        }
        else
            noImpl();
    }

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> sinh(const Diff<ScalarType, Mode, Order>& s) {
        using ResultType = Diff<ScalarType, Mode, Order>;
        if constexpr (Mode == DiffMode::Forward) {
            using GradType = typename ResultType::GradType;
            return ResultType(sinh(s.getValue()), s.getGrad() * cosh(GradType(s)));
        }
        else
            noImpl();
    }

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> tanh(const Diff<ScalarType, Mode, Order>& s) {
        using ResultType = Diff<ScalarType, Mode, Order>;
        if constexpr (Mode == DiffMode::Forward) {
            using GradType = typename ResultType::GradType;
            const GradType v = tanh(GradType(s));
            return ResultType(v.getValue(), (ScalarType(1) - square(v)) * s.getGrad());
        }
        else {
            auto& tracer = DiffTracer<ScalarType, Order>::getInstance();
            const auto result = tracer.pushOperation(tanh(s.getValue()), ExprType::Tanh);
            tracer.pushOperand(s);
            return result;
        }
    }

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> lncosh(const Diff<ScalarType, Mode, Order>& s) noexcept {
        using ResultType = Diff<ScalarType, Mode, Order>;
        if constexpr (Mode == DiffMode::Forward) {
            using GradType = typename ResultType::GradType;
            return ResultType(lncosh(s.getValue()), s.getGrad() * tanh(GradType(s)));
        }
        else
            noImpl();
    }
}
