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
    template<Scalar T, size_t Size>
    [[nodiscard]] __host__ __device__ inline SIMD<T, Size> abs(const SIMD<T, Size>& x) {
        if constexpr (T::Option == Float16)
            return SIMD<T, Size>(__habs2(x.toMachine()));
        else if constexpr (IsHost())
            return SIMD<T, Size>(abs(x.toMachine()));
        else
            noImpl();
    }

    template<Scalar T, size_t Size>
    [[nodiscard]] inline SIMD<T, Size> square(const SIMD<T, Size>& x) {
        return x * x;
    }

    template<Scalar T, size_t Size>
    [[nodiscard]] inline auto reciprocal(const SIMD<T, Size>& x) {
        using ResultType = SIMD<T, Size>;
        return ResultType(1) / x;
    }

    template<Scalar T, size_t Size>
    [[nodiscard]] inline SIMD<T, Size> sqrt(const SIMD<T, Size>& x) {
        if constexpr (!T::isDifferentiable)
            return SIMD<T, Size>(sqrt(x.toMachine()));
        else {
            using PlainSIMD = SIMD<typename T::ValueType, Size>;
            using TracerType = typename T::TracerType;
            auto& tracer = TracerType::getInstance();
            const PlainSIMD values(sqrt(x.toMachine()));
            const size_t newHeadNode = tracer.pushOperation(values, ExprType::Sqrt);
            for (size_t i = 0; i < Size; ++i)
                tracer.pushOperand(T(x.value_ptr() + i, x.grad_ptr() + i));
            return {values, newHeadNode};
        }
    }

    template<Scalar T, size_t Size>
    [[nodiscard]] inline SIMD<T, Size> cbrt(const SIMD<T, Size>& x) {
        if constexpr (!T::isDifferentiable)
            return SIMD<T, Size>(cbrt(x.toMachine()));
        else {
            using PlainSIMD = SIMD<typename T::ValueType, Size>;
            using TracerType = typename T::TracerType;
            auto& tracer = TracerType::getInstance();
            const PlainSIMD values(cbrt(x.toMachine()));
            const size_t newHeadNode = tracer.pushOperation(values, ExprType::Cbrt);
            for (size_t i = 0; i < Size; ++i)
                tracer.pushOperand(T(x.value_ptr() + i, x.grad_ptr() + i));
            return {values, newHeadNode};
        }
    }

    template<Scalar T, size_t Size>
    [[nodiscard]] inline auto ln(const SIMD<T, Size>& x) {
        return SIMD<T, Size>(log(x.toMachine()));
    }

    template<Scalar T, size_t Size>
    [[nodiscard]] inline auto ln1p(const SIMD<T, Size>& x) {
        return SIMD<T, Size>(log1p(x.toMachine()));
    }

    template<Scalar T, size_t Size>
    [[nodiscard]] inline auto exp(const SIMD<T, Size>& x) {
        return SIMD<T, Size>(exp(x.toMachine()));
    }

    template<Scalar T, size_t Size>
    inline static void sincos(const SIMD<T, Size>& x, SIMD<T, Size>& s, SIMD<T, Size>& c) {
        sincos(x.toMachine(), s.toMachine(), c.toMachine());
        if constexpr (T::isDifferentiable) {
            using TracerType = typename T::TracerType;
            auto& tracer = TracerType::getInstance();
            const T sinHeadNode = tracer.pushOperation(s, ExprType::Sin);
            for (size_t i = 0; i < Size; ++i)
                tracer.pushOperand(T(x.value_ptr() + i, x.grad_ptr() + i));
            s = SIMD<T, Size>(s, sinHeadNode);

            const T cosHeadNode = tracer.pushOperation(c, ExprType::Cos);
            for (size_t i = 0; i < Size; ++i)
                tracer.pushOperand(T(x.value_ptr() + i, x.grad_ptr() + i));
            c = SIMD<T, Size>(c, cosHeadNode);
        }
    }

    template<Scalar T, size_t Size>
    [[nodiscard]] inline SIMD<T, Size> lncosh(const SIMD<T, Size>& x) noexcept {
        using Pack = SIMD<T, Size>;
        const auto x1 = abs(x);
        return x1 + ln1p(exp(-x1 * Pack(2))) - Pack(M_LN2);
    }
}
