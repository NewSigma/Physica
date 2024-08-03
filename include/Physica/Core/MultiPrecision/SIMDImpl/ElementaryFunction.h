/*
 * Copyright 2023 Weibo He.
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
            using PlainSIMD = SIMD<typename ScalarType::PlainScalar, Size>;
            using TracerType = typename ScalarType::TracerType;
            auto& tracer = TracerType::getInstance();
            const PlainSIMD values(abs(x.getImpl()));
            const ScalarType newHeadNode = tracer.pushOperation(values, ExpressionType::Abs);
            for (size_t i = 0; i < Size; ++i)
                tracer.pushOperand(ScalarType(x.value_ptr() + i, x.grad_ptr() + i));
            return {values, newHeadNode};
        }
    }

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> square(const SIMD<ScalarType, Size>& x) {
        if constexpr (!ScalarType::isDifferentiable)
            return x * x;
        else {
            using PlainSIMD = SIMD<typename ScalarType::PlainScalar, Size>;
            using TracerType = typename ScalarType::TracerType;
            auto& tracer = TracerType::getInstance();
            const PlainSIMD values(square(x.getImpl()));
            const size_t newHeadNode = tracer.pushOperation(values, ExpressionType::Square);
            for (size_t i = 0; i < Size; ++i)
                tracer.pushOperand(ScalarType(x.value_ptr() + i, x.grad_ptr() + i));
            return {values, newHeadNode};
        }
    }

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> sqrt(const SIMD<ScalarType, Size>& x) {
        if constexpr (!ScalarType::isDifferentiable)
            return SIMD<ScalarType, Size>(sqrt(x.getImpl()));
        else {
            using PlainSIMD = SIMD<typename ScalarType::PlainScalar, Size>;
            using TracerType = typename ScalarType::TracerType;
            auto& tracer = TracerType::getInstance();
            const PlainSIMD values(sqrt(x.getImpl()));
            const size_t newHeadNode = tracer.pushOperation(values, ExpressionType::Sqrt);
            for (size_t i = 0; i < Size; ++i)
                tracer.pushOperand(ScalarType(x.value_ptr() + i, x.grad_ptr() + i));
            return {values, newHeadNode};
        }
    }

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> cbrt(const SIMD<ScalarType, Size>& x) {
        if constexpr (!ScalarType::isDifferentiable)
            return SIMD<ScalarType, Size>(cbrt(x.getImpl()));
        else {
            using PlainSIMD = SIMD<typename ScalarType::PlainScalar, Size>;
            using TracerType = typename ScalarType::TracerType;
            auto& tracer = TracerType::getInstance();
            const PlainSIMD values(cbrt(x.getImpl()));
            const size_t newHeadNode = tracer.pushOperation(values, ExpressionType::Cbrt);
            for (size_t i = 0; i < Size; ++i)
                tracer.pushOperand(ScalarType(x.value_ptr() + i, x.grad_ptr() + i));
            return {values, newHeadNode};
        }
    }

    template<class ScalarType, size_t Size>
    inline static void sincos(const SIMD<ScalarType, Size>& x, SIMD<ScalarType, Size>& s, SIMD<ScalarType, Size>& c) {
        sincos(x.getImpl(), s.getImpl(), c.getImpl());
        if constexpr (ScalarType::isDifferentiable) {
            using TracerType = typename ScalarType::TracerType;
            auto& tracer = TracerType::getInstance();
            const ScalarType sinHeadNode = tracer.pushOperation(s, ExpressionType::Sin);
            for (size_t i = 0; i < Size; ++i)
                tracer.pushOperand(ScalarType(x.value_ptr() + i, x.grad_ptr() + i));
            s = SIMD<ScalarType, Size>(s, sinHeadNode);

            const ScalarType cosHeadNode = tracer.pushOperation(c, ExpressionType::Cos);
            for (size_t i = 0; i < Size; ++i)
                tracer.pushOperand(ScalarType(x.value_ptr() + i, x.grad_ptr() + i));
            c = SIMD<ScalarType, Size>(c, cosHeadNode);
        }
    }
}
