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

#include "Math.h"

namespace Physica::Core {
    template<Scalar T>
    __host__ __device__ inline CoDiff<T> abs(T&& x) requires(Diffable<T>) {
        using ScalarType = std::remove_reference_t<T>::ScalarType;
        if constexpr (ForwardDiff<T>)
            co_return ScalarType(abs(x.value()), x.value().isPositive() ? x.grad() : -x.grad());
        else {
            LazyReverse<T&&> x_ = std::forward<T>(x);
            auto result = co_yield abs(x_.value());
            auto& g = result.grad();
            if (!g.isZero())
                x_.reverse(x_.isPositive() ? g : -g);
        }
    }

    template<Scalar T, int Order>
    __host__ __device__ inline auto abs(const ScalarRef<Diff<T, DiffMode::Forward, Order>>& x) {
        return abs(Diff<T, DiffMode::Forward, Order>(x));
    }

    template<Scalar T>
    inline CoDiff<T> relu(T&& x) requires(Diffable<T>) {
        using ScalarType = std::remove_reference_t<T>::ScalarType;
        using ValueType = ScalarType::ValueType;
        if constexpr (ForwardDiff<T>) {
            co_return ScalarType(relu(x.value()), x.value().isPositive() ? x.grad() : T(0));
        }
        else {
            LazyReverse<T&&> x_ = std::forward<T>(x);
            auto result = co_yield relu(x_.value());
            auto& g = result.grad();
            if (!g.isZero())
                x_.reverse(x_.isPositive() ? g : ValueType(0));
        }
    }

    template<Scalar T>
    __host__ __device__ inline CoDiff<T> square(T&& x) requires(Diffable<T>) {
        using ScalarType = std::remove_reference_t<T>::ScalarType;
        using ValueType = ScalarType::ValueType;
        if constexpr (ForwardDiff<T>) {
            using GradType = ScalarType::GradType;
            co_return ScalarType(square(x.value()), GradType(ValueType(2) * x * x.grad()));
        }
        else {
            LazyReverse<T&&> x_ = std::forward<T>(x);
            auto result = co_yield square(x_.value());
            auto& g = result.grad();
            if (!g.isZero())
                x_.reverse(ValueType(2) * x_.value() * g);
        }
    }

    template<Scalar T>
    __host__ __device__ inline CoDiff<T> reciprocal(T&& x) requires(Diffable<T>) {
        using ScalarType = std::remove_reference_t<T>::ScalarType;
        if constexpr (ForwardDiff<T>) {
            using GradType = ScalarType::GradType;
            const GradType v = reciprocal(GradType(x));
            co_return ScalarType(v.value(), -x.grad() * square(v));
        }
        else {
            LazyReverse<T&&> x_ = std::forward<T>(x);
            auto result = co_yield reciprocal(x_.value());
            auto& g = result.grad();
            if (!g.isZero())
                x_.reverse(-square(result.value()) * g);
        }
    }

    template<Scalar T>
    __host__ __device__ CoDiff<T> sqrt(T&& x) requires(Diffable<T>) {
        using ScalarType = std::remove_reference_t<T>::ScalarType;
        using ValueType = ScalarType::ValueType;
        if constexpr (ForwardDiff<T>) {
            using GradType = ScalarType::GradType;
            const GradType v = sqrt(GradType(x));
            co_return ScalarType(v.value(), ValueType(0.5) * x.grad() / v);
        }
        else {
            LazyReverse<T&&> x_ = std::forward<T>(x);
            auto result = co_yield sqrt(x_.value());
            auto& g = result.grad();
            if (!g.isZero())
                x_.reverse(reciprocal(result.value()) * ValueType(0.5));
        }
    }

    template<Scalar T>
    CoDiff<T> cbrt(T&& x) requires(Diffable<T>) {
        using ScalarType = std::remove_reference_t<T>::ScalarType;
        using ValueType = ScalarType::ValueType;
        if constexpr (ForwardDiff<T>) {
            using GradType = ScalarType::GradType;
            const GradType v = cbrt(GradType(x));
            co_return ScalarType(v.value(), ValueType(1.0 / 3) * v * x.grad() / GradType(x));
        }
        else {
            LazyReverse<T&&> x_ = std::forward<T>(x);
            auto result = co_yield cbrt(x_.value());
            auto& g = result.grad();
            if (!g.isZero()) {
                const auto x2_3 = result.value() / x_.value();
                x_.reverse((ValueType(1.0 / 3) * g) * x2_3);
            }
        }
    }

    template<Scalar T>
    CoDiff<T> ln(T&& x) requires(Diffable<T>) {
        using ScalarType = std::remove_reference_t<T>::ScalarType;
        if constexpr (ForwardDiff<T>) {
            using GradType = ScalarType::GradType;
            co_return ScalarType(ln(x.value()), x.grad() / GradType(x));
        }
        else {
            LazyReverse<T&&> x_ = std::forward<T>(x);
            auto result = co_yield ln(x_.value());
            auto& g = result.grad();
            if (!g.isZero())
                x_.reverse(result.value() * g / x_.value());
        }
    }

    template<Scalar T, int Order>
    auto ln1p(const Diff<T, DiffMode::Forward, Order>& x) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        using GradType = ResultType::GradType;
        return ResultType(ln1p(x.value()), x.grad() / (T(1) + GradType(x)));
    }

    template<Scalar T, int Order>
    auto ln1pexp(const Diff<T, DiffMode::Forward, Order>& x) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        using GradType = ResultType::GradType;
        const GradType x1 = GradType(x);
        return ResultType(ln1pexp(x.value()), reciprocal(GradType(1) + exp(-x1)) * x.grad());
    }

    template<Scalar T>
    CoDiff<T> exp(T&& x) requires(Diffable<T>) {
        using ScalarType = std::remove_reference_t<T>::ScalarType;
        if constexpr (ForwardDiff<T>) {
            using GradType = ScalarType::GradType;
            const GradType v = exp(GradType(x));
            co_return ScalarType(v.value(), v * x.grad());
        }
        else {
            LazyReverse<T&&> x_ = std::forward<T>(x);
            auto result = co_yield exp(x_.value());
            auto& g = result.grad();
            if (!g.isZero())
                x_.reverse(result.value() * g);
        }
    }

    template<Scalar T, int Order>
    auto pow(const Diff<T, DiffMode::Forward, Order>& x, const T& a) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        using GradType = ResultType::GradType;
        constexpr int GradOrder = GradType::Order;
        const auto y = pow(GradType(x), a);
        return ResultType(y.value(), x.grad() * y / x.template mask<GradOrder>() * a);
    }

    template<Scalar T>
    CoDiff<T> cos(T&& x) requires(Diffable<T>) {
        using ScalarType = std::remove_reference_t<T>::ScalarType;
        if constexpr (ForwardDiff<T>) {
            using GradType = ScalarType::GradType;
            GradType sin_value, cos_value;
            sincos(GradType(x), sin_value, cos_value);
            co_return ScalarType(cos_value.value(), -sin_value * x.grad());
        }
        else {
            using ValueType = ScalarType::ValueType;
            LazyReverse<T&&> x_ = std::forward<T>(x);
            ValueType c, s;
            sincos(x_.value(), s, c);
            auto result = co_yield c;
            auto& g = result.grad();
            if (!g.isZero())
                x_.reverse(-s * g);
        }
    }

    template<Scalar T>
    CoDiff<T> sin(T&& x) requires(Diffable<T>) {
        using ScalarType = std::remove_reference_t<T>::ScalarType;
        if constexpr (ForwardDiff<T>) {
            using GradType = ScalarType::GradType;
            GradType sin_value, cos_value;
            sincos(GradType(x), sin_value, cos_value);
            co_return ScalarType(sin_value.value(), cos_value * x.grad());
        }
        else {
            using ValueType = ScalarType::ValueType;
            LazyReverse<T&&> x_ = std::forward<T>(x);
            ValueType c, s;
            sincos(x_.value(), s, c);
            auto result = co_yield s;
            auto& g = result.grad();
            if (!g.isZero())
                x_.reverse(c * g);
        }
    }

    template<Scalar T, Scalar U>
    CoDiff<void> sincos(const T& x, U&& sin_result, U&& cos_result) {
        if constexpr (ForwardDiff<T>) {
            using GradType = T::GradType;
            GradType sin_value, cos_value;
            sincos(GradType(x), sin_value, cos_value);
            sin_result = T(sin_value, cos_value * x.grad());
            cos_result = T(cos_value, -sin_value * x.grad());
            co_return;
        }
        else {
            using ValueType = T::ValueType;
            ValueType s, c;
            sincos(x.value(), s, c);
            LazyReverse<U&&> sin_ = std::forward<U>(sin_result);
            LazyReverse<U&&> cos_ = std::forward<U>(cos_result);
            sin_ = s;
            cos_ = c;
            co_await std::suspend_always{};
        }
    }

    template<Scalar T, int Order>
    auto tan(const Diff<T, DiffMode::Forward, Order>& x) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        using GradType = ResultType::GradType;
        return ResultType(tan(x.value()), x.grad() * square(sec(GradType(x))));
    }

    template<Scalar T, int Order>
    auto sec(const Diff<T, DiffMode::Forward, Order>& x) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        using GradType = ResultType::GradType;
        const auto x1 = GradType(x);
        const auto v = sec(x1);
        return ResultType(v.value(), x.grad() * v * tan(x1));
    }

    template<Scalar T, int Order>
    auto csc(const Diff<T, DiffMode::Forward, Order>& x) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        using GradType = ResultType::GradType;
        const auto x1 = GradType(x);
        const auto v = csc(x1);
        return ResultType(v.value(), -x.grad() * v * cot(x1));
    }

    template<Scalar T, int Order>
    auto cot(const Diff<T, DiffMode::Forward, Order>& x) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        using GradType = ResultType::GradType;
        return ResultType(cot(x.value()), -x.grad() * square(csc(GradType(x))));
    }

    template<Scalar T, int Order>
    auto arccos(const Diff<T, DiffMode::Forward, Order>& x) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        using GradType = ResultType::GradType;
        return ResultType(arccos(x.value()), -x.grad() / sqrt(T(1) - square(GradType(x))));
    }

    template<Scalar T, int Order>
    auto cosh(const Diff<T, DiffMode::Forward, Order>& x) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        using GradType = ResultType::GradType;
        return ResultType(cosh(x.value()), x.grad() * sinh(GradType(x)));
    }

    template<Scalar T, int Order>
    auto sinh(const Diff<T, DiffMode::Forward, Order>& x) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        using GradType = ResultType::GradType;
        return ResultType(sinh(x.value()), x.grad() * cosh(GradType(x)));
    }

    template<Scalar T, int Order>
    auto tanh(const Diff<T, DiffMode::Forward, Order>& x) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        using GradType = ResultType::GradType;
        const GradType v = tanh(GradType(x));
        return ResultType(v.value(), (T(1) - square(v)) * x.grad());
    }

    template<Scalar T, int Order>
    auto lncosh(const Diff<T, DiffMode::Forward, Order>& x) noexcept {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        using GradType = ResultType::GradType;
        return ResultType(lncosh(x.value()), x.grad() * tanh(GradType(x)));
    }
}
