/*
 * Copyright 2024 Weibo He.
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
    template<Scalar T, DiffMode Mode, int Order, size_t Size>
    [[nodiscard]] inline auto abs(const SIMD<Diff<T, Mode, Order>, Size>& x) {
        using ResultType = SIMD<Diff<T, Mode, Order>, Size>;
        if constexpr (Mode == DiffMode::Forward) {
            using GradPacket = ResultType::GradType;
            return ResultType(abs(x.value()), GradPacket::select(x.value().isPositive(), x.grad(), -x.grad()));
        }
        else {
            using PlainSIMD = SIMD<T, Size>;
            using TracerType = T::TracerType;
            auto& tracer = TracerType::getInstance();
            const PlainSIMD values(abs(x.toMachine()));
            const auto newHeadNode = tracer.pushOperation(values, ExprType::Abs);
            for (size_t i = 0; i < Size; ++i)
                tracer.pushOperand(T(x.value_ptr() + i, x.grad_ptr() + i));
            return ResultType(values, newHeadNode);
        }
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Size>
    [[nodiscard]] inline auto square(const SIMD<Diff<T, Mode, Order>, Size>& x) {
        using ResultType = SIMD<Diff<T, Mode, Order>, Size>;
        if constexpr (Mode == DiffMode::Forward) {
            using GradPacket = ResultType::GradType;
            return ResultType(square(x.value()), GradPacket(x) * x.grad() * T(2));
        }
        else {
            using PlainSIMD = SIMD<T, Size>;
            using TracerType = T::TracerType;
            auto& tracer = TracerType::getInstance();
            const PlainSIMD values(square(x.toMachine()));
            const auto newHeadNode = tracer.pushOperation(values, ExprType::Square);
            for (size_t i = 0; i < Size; ++i)
                tracer.pushOperand(T(x.value_ptr() + i, x.grad_ptr() + i));
            return ResultType(values, newHeadNode);
        }
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Size>
    [[nodiscard]] inline auto reciprocal(const SIMD<Diff<T, Mode, Order>, Size>& x) {
        static_assert(Mode != DiffMode::Reverse, "[Error]: Not implemented");
        using ResultType = SIMD<Diff<T, Mode, Order>, Size>;
        using GradPacket = ResultType::GradType;
        const auto y = reciprocal(GradPacket(x));
        return ResultType(y.value(), -x.grad() * square(y));
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Size>
    [[nodiscard]] inline auto ln(const SIMD<Diff<T, Mode, Order>, Size>& x) {
        static_assert(Mode != DiffMode::Reverse, "[Error]: Not implemented");
        using ResultType = SIMD<Diff<T, Mode, Order>, Size>;
        using GradPacket = ResultType::GradType;
        return ResultType(ln(x.value()), reciprocal(GradPacket(x)) * x.grad());
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Size>
    [[nodiscard]] inline auto ln1p(const SIMD<Diff<T, Mode, Order>, Size>& x) {
        static_assert(Mode != DiffMode::Reverse, "[Error]: Not implemented");
        using ResultType = SIMD<Diff<T, Mode, Order>, Size>;
        using GradPacket = ResultType::GradType;
        return ResultType(ln1p(x.value()), reciprocal(GradPacket(1) + GradPacket(x)) * x.grad());
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Size>
    [[nodiscard]] inline auto exp(const SIMD<Diff<T, Mode, Order>, Size>& x) {
        static_assert(Mode != DiffMode::Reverse, "[Error]: Not implemented");
        using ResultType = SIMD<Diff<T, Mode, Order>, Size>;
        using GradPacket = ResultType::GradType;
        const auto y = exp(GradPacket(x));
        return ResultType(y.value(), y * x.grad());
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Size>
    [[nodiscard]] inline auto sincos(
            const SIMD<Diff<T, Mode, Order>, Size>& x,
            SIMD<Diff<T, Mode, Order>, Size>& s,
            SIMD<Diff<T, Mode, Order>, Size>& c) {
        static_assert(Mode == DiffMode::Reverse, "[Error]: Not implemented");
        using ScalarType = Diff<T, Mode, Order>;
        Physica::sincos(x.toMachine(), s.toMachine(), c.toMachine());
        if constexpr (Mode == DiffMode::Reverse) {
            using TracerType = ScalarType::TracerType;
            auto& tracer = TracerType::getInstance();
            const auto sinHeadNode = tracer.pushOperation(s, ExprType::Sin);
            for (size_t i = 0; i < Size; ++i)
                tracer.pushOperand(ScalarType(x.value_ptr() + i, x.grad_ptr() + i));
            s = SIMD<ScalarType, Size>(s, sinHeadNode);

            const auto cosHeadNode = tracer.pushOperation(c, ExprType::Cos);
            for (size_t i = 0; i < Size; ++i)
                tracer.pushOperand(ScalarType(x.value_ptr() + i, x.grad_ptr() + i));
            c = SIMD<ScalarType, Size>(c, cosHeadNode);
        }
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Size>
    [[nodiscard]] inline auto tanh(const SIMD<Diff<T, Mode, Order>, Size>& x) {
        static_assert(Mode != DiffMode::Reverse, "[Error]: Not implemented");
        using ResultType = SIMD<Diff<T, Mode, Order>, Size>;
        using GradPacket = ResultType::GradType;
        const GradPacket y = tanh(GradPacket(x));
        return ResultType(y.value(), (GradPacket(1) - square(y)) * x.grad());
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Size>
    [[nodiscard]] inline auto lncosh(const SIMD<Diff<T, Mode, Order>, Size>& x) {
        static_assert(Mode != DiffMode::Reverse, "[Error]: Not implemented");
        using ResultType = SIMD<Diff<T, Mode, Order>, Size>;
        using GradPacket = ResultType::GradType;
        return ResultType(lncosh(x.value()), tanh(GradPacket(x)) * x.grad());
    }
}
