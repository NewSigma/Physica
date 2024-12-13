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
    template<Scalar T, int Order>
    __host__ __device__ inline auto abs(const Diff<T, DiffMode::Forward, Order>& x) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        return ResultType(abs(x.getValue()), x.getValue().isPositive() ? x.getGrad() : -x.getGrad());
    }

    template<Scalar T, int Order>
    __host__ __device__ inline auto abs(const ScalarRef<Diff<T, DiffMode::Forward, Order>>& x) {
        return abs(Diff<T, DiffMode::Forward, Order>(x));
    }

    template<Scalar T, int Order>
    inline auto relu(const Diff<T, DiffMode::Forward, Order>& x) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        return ResultType(relu(x.getValue()), x.getValue().isPositive() ? x.getGrad() : T(0));
    }

    template<Scalar T>
    __host__ __device__ inline CoDiff<T> square(T&& x) requires(Diffable<T>) {
        using ScalarType = std::remove_reference_t<T>::ScalarType;
        using ValueType = ScalarType::ValueType;
        if constexpr (ForwardDiff<T>) {
            using GradType = ScalarType::GradType;
            co_return ScalarType(square(x.getValue()), GradType(ValueType(2) * x * x.getGrad()));
        }
        else {
            LazyReverse<T&&> x_ = std::forward<T>(x);
            auto result = co_yield square(x_.getValue());
            auto& g = result.getGrad();
            if (!g.isZero())
                x_.reverse(ValueType(2) * x_.getValue() * g);
        }
    }

    template<Scalar T>
    __host__ __device__ inline CoDiff<T> reciprocal(T&& x) requires(Diffable<T>) {
        using ScalarType = std::remove_reference_t<T>::ScalarType;
        if constexpr (ForwardDiff<T>) {
            using GradType = ScalarType::GradType;
            const GradType v = reciprocal(GradType(x));
            co_return ScalarType(v.getValue(), -x.getGrad() * square(v));
        }
        else {
            LazyReverse<T&&> x_ = std::forward<T>(x);
            auto result = co_yield reciprocal(x_.getValue());
            auto& g = result.getGrad();
            if (!g.isZero())
                x_.reverse(-square(result.getValue()) * g);
        }
    }

    template<Scalar T>
    __host__ __device__ CoDiff<T> sqrt(T&& x) requires(Diffable<T>) {
        using ScalarType = std::remove_reference_t<T>::ScalarType;
        using ValueType = ScalarType::ValueType;
        if constexpr (ForwardDiff<T>) {
            using GradType = ScalarType::GradType;
            const GradType v = sqrt(GradType(x));
            co_return ScalarType(v.getValue(), ValueType(0.5) * x.getGrad() / v);
        }
        else {
            LazyReverse<T&&> x_ = std::forward<T>(x);
            auto result = co_yield sqrt(x_.getValue());
            auto& g = result.getGrad();
            if (!g.isZero())
                x_.reverse(reciprocal(result.getValue()) * ValueType(0.5));
        }
    }

    template<Scalar T>
    CoDiff<T> cbrt(T&& x) requires(Diffable<T>) {
        using ScalarType = std::remove_reference_t<T>::ScalarType;
        using ValueType = ScalarType::ValueType;
        if constexpr (ForwardDiff<T>) {
            using GradType = ScalarType::GradType;
            const GradType v = cbrt(GradType(x));
            co_return ScalarType(v.getValue(), ValueType(1.0 / 3) * v * x.getGrad() / GradType(x));
        }
        else {
            LazyReverse<T&&> x_ = std::forward<T>(x);
            auto result = co_yield cbrt(x_.getValue());
            auto& g = result.getGrad();
            if (!g.isZero()) {
                const auto x2_3 = result.getValue() / x_.getValue();
                x_.reverse((ValueType(1.0 / 3) * g) * x2_3);
            }
        }
    }

    template<Scalar T, int Order>
    auto ln(const Diff<T, DiffMode::Forward, Order>& x) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        using GradType = ResultType::GradType;
        return ResultType(ln(x.getValue()), x.getGrad() / GradType(x));
    }

    template<Scalar T, int Order>
    auto ln1p(const Diff<T, DiffMode::Forward, Order>& x) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        using GradType = ResultType::GradType;
        return ResultType(ln1p(x.getValue()), x.getGrad() / (T(1) + GradType(x)));
    }

    template<Scalar T>
    CoDiff<T> exp(T&& x) requires(Diffable<T>) {
        using ScalarType = std::remove_reference_t<T>::ScalarType;
        if constexpr (ForwardDiff<T>) {
            using GradType = ScalarType::GradType;
            const GradType v = exp(GradType(x));
            co_return ScalarType(v.getValue(), v * x.getGrad());
        }
        else {
            LazyReverse<T&&> x_ = std::forward<T>(x);
            auto result = co_yield exp(x_.getValue());
            auto& g = result.getGrad();
            if (!g.isZero())
                x_.reverse(result.getValue() * g);
        }
    }

    template<Scalar T, int Order>
    auto pow(const Diff<T, DiffMode::Forward, Order>& x, const T& a) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        using GradType = ResultType::GradType;
        constexpr int GradOrder = GradType::Order;
        const auto y = pow(GradType(x), a);
        return ResultType(y.getValue(), x.getGrad() * y / x.template mask<GradOrder>() * a);
    }

    template<Scalar T>
    CoDiff<T> cos(T&& x) requires(Diffable<T>) {
        using ScalarType = std::remove_reference_t<T>::ScalarType;
        if constexpr (ForwardDiff<T>) {
            using GradType = ScalarType::GradType;
            GradType sin_value, cos_value;
            sincos(GradType(x), sin_value, cos_value);
            co_return ScalarType(cos_value.getValue(), -sin_value * x.getGrad());
        }
        else {
            using ValueType = ScalarType::ValueType;
            LazyReverse<T&&> x_ = std::forward<T>(x);
            ValueType c, s;
            sincos(x_.getValue(), s, c);
            auto result = co_yield c;
            auto& g = result.getGrad();
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
            co_return ScalarType(sin_value.getValue(), cos_value * x.getGrad());
        }
        else {
            using ValueType = ScalarType::ValueType;
            LazyReverse<T&&> x_ = std::forward<T>(x);
            ValueType c, s;
            sincos(x_.getValue(), s, c);
            auto result = co_yield s;
            auto& g = result.getGrad();
            if (!g.isZero())
                x_.reverse(c * g);
        }
    }

    template<Scalar T, int Order>
    void sincos(const Diff<T, DiffMode::Forward, Order>& x, Diff<T, DiffMode::Forward, Order>& sin_result, Diff<T, DiffMode::Forward, Order>& cos_result) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        using GradType = ResultType::GradType;
        GradType sin_value, cos_value;
        sincos(GradType(x), sin_value, cos_value);
        sin_result = ResultType(sin_value, cos_value * x.getGrad());
        cos_result = ResultType(cos_value, -sin_value * x.getGrad());
    }

    template<Scalar T, int Order>
    auto tan(const Diff<T, DiffMode::Forward, Order>& x) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        using GradType = ResultType::GradType;
        return ResultType(tan(x.getValue()), x.getGrad() * square(sec(GradType(x))));
    }

    template<Scalar T, int Order>
    auto sec(const Diff<T, DiffMode::Forward, Order>& x) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        using GradType = ResultType::GradType;
        const auto x1 = GradType(x);
        const auto v = sec(x1);
        return ResultType(v.getValue(), x.getGrad() * v * tan(x1));
    }

    template<Scalar T, int Order>
    auto csc(const Diff<T, DiffMode::Forward, Order>& x) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        using GradType = ResultType::GradType;
        const auto x1 = GradType(x);
        const auto v = csc(x1);
        return ResultType(v.getValue(), -x.getGrad() * v * cot(x1));
    }

    template<Scalar T, int Order>
    auto cot(const Diff<T, DiffMode::Forward, Order>& x) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        using GradType = ResultType::GradType;
        return ResultType(cot(x.getValue()), -x.getGrad() * square(csc(GradType(x))));
    }

    template<Scalar T, int Order>
    auto arccos(const Diff<T, DiffMode::Forward, Order>& x) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        using GradType = ResultType::GradType;
        return ResultType(arccos(x.getValue()), -x.getGrad() / sqrt(T(1) - square(GradType(x))));
    }

    template<Scalar T, int Order>
    auto cosh(const Diff<T, DiffMode::Forward, Order>& x) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        using GradType = ResultType::GradType;
        return ResultType(cosh(x.getValue()), x.getGrad() * sinh(GradType(x)));
    }

    template<Scalar T, int Order>
    auto sinh(const Diff<T, DiffMode::Forward, Order>& x) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        using GradType = ResultType::GradType;
        return ResultType(sinh(x.getValue()), x.getGrad() * cosh(GradType(x)));
    }

    template<Scalar T, int Order>
    auto tanh(const Diff<T, DiffMode::Forward, Order>& x) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        using GradType = ResultType::GradType;
        const GradType v = tanh(GradType(x));
        return ResultType(v.getValue(), (T(1) - square(v)) * x.getGrad());
    }

    template<Scalar T, int Order>
    auto lncosh(const Diff<T, DiffMode::Forward, Order>& x) noexcept {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        using GradType = ResultType::GradType;
        return ResultType(lncosh(x.getValue()), x.getGrad() * tanh(GradType(x)));
    }
}
