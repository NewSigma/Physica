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
    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> abs(const SIMD<ScalarType, Size>& x) {
        if constexpr (!ScalarType::isDifferentiable)
            return SIMD<ScalarType, Size>(abs(x.getImpl()));
        else {
            using PlainScalar = typename ScalarType::PlainScalar;
            auto& tracer = DiffTracer<PlainScalar>::getInstance();
            const auto temp = abs(x.getImpl());
            const size_t headIndex = x.getHeadTraceIndex();
            const size_t newHeadIndex = tracer.pushOperation(temp, ExpressionType::Abs);
            for (size_t i = 0; i < Size; ++i)
                tracer.pushOperand(headIndex + i);
            return {temp, newHeadIndex};
        }
    }

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> square(const SIMD<ScalarType, Size>& x) {
        if constexpr (!ScalarType::isDifferentiable)
            return x * x;
        else {
            using PlainScalar = typename ScalarType::PlainScalar;
            auto& tracer = DiffTracer<PlainScalar>::getInstance();
            const auto temp = square(x.getImpl());
            const size_t headIndex = x.getHeadTraceIndex();
            const size_t newHeadIndex = tracer.pushOperation(temp, ExpressionType::Square);
            for (size_t i = 0; i < Size; ++i)
                tracer.pushOperand(headIndex + i);
            return {temp, newHeadIndex};
        }
    }

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> sqrt(const SIMD<ScalarType, Size>& x) {
        if constexpr (!ScalarType::isDifferentiable)
            return SIMD<ScalarType, Size>(sqrt(x.getImpl()));
        else {
            using PlainScalar = typename ScalarType::PlainScalar;
            auto& tracer = DiffTracer<PlainScalar>::getInstance();
            const auto temp = sqrt(x.getImpl());
            const size_t headIndex = x.getHeadTraceIndex();
            const size_t newHeadIndex = tracer.pushOperation(temp, ExpressionType::Sqrt);
            for (size_t i = 0; i < Size; ++i)
                tracer.pushOperand(headIndex + i);
            return {temp, newHeadIndex};
        }
    }

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> cbrt(const SIMD<ScalarType, Size>& x) {
        if constexpr (!ScalarType::isDifferentiable)
            return SIMD<ScalarType, Size>(cbrt(x.getImpl()));
        else {
            using PlainScalar = typename ScalarType::PlainScalar;
            auto& tracer = DiffTracer<PlainScalar>::getInstance();
            const auto temp = cbrt(x.getImpl());
            const size_t headIndex = x.getHeadTraceIndex();
            const size_t newHeadIndex = tracer.pushOperation(temp, ExpressionType::Cbrt);
            for (size_t i = 0; i < Size; ++i)
                tracer.pushOperand(headIndex + i);
            return {temp, newHeadIndex};
        }
    }

    template<class ScalarType, size_t Size>
    inline static void sincos(const SIMD<ScalarType, Size>& x, SIMD<ScalarType, Size>& s, SIMD<ScalarType, Size>& c) {
        sincos(x.getImpl(), s.getImpl(), c.getImpl());
        if constexpr (ScalarType::isDifferentiable) {
            using PlainScalar = typename ScalarType::PlainScalar;
            auto& tracer = DiffTracer<PlainScalar>::getInstance();
            const size_t headIndex = x.getHeadTraceIndex();
            const ScalarType sinHeadTrace = tracer.pushOperation(s, ExpressionType::Sin);
            for (size_t i = 0; i < Size; ++i)
                tracer.pushOperand(headIndex + i);
            s = SIMD<ScalarType, Size>(s, sinHeadTrace);

            const ScalarType cosHeadTrace = tracer.pushOperation(c, ExpressionType::Cos);
            for (size_t i = 0; i < Size; ++i)
                tracer.pushOperand(headIndex + i);
            c = SIMD<ScalarType, Size>(c, cosHeadTrace);
        }
    }
}
