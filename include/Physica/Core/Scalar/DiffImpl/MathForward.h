/*
 * Copyright 2025 Weibo He.
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

#include "Physica/Core/Scalar/Diff.h"

namespace Physica {
    template<Scalar T>
    __host__ __device__ auto abs(const T& x) requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        return ResultType(abs(x.value()), x.value().isPositive() ? x.grad() : -x.grad());
    }

    template<Scalar T>
    __host__ __device__ auto relu(const T& x) requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        using GradType = T::GradType;
        return ResultType(relu(x.value()), x.value().isPositive() ? x.grad() : GradType(0));
    }

    template<Scalar T>
    __host__ __device__ auto square(const T& x) requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        using Tv = T::ValueType;
        constexpr int GradOrder = T::GradType::Order;
        return ResultType(square(x.value()), Tv(2) * x.template mask<GradOrder>() * x.grad());
    }

    template<Scalar T>
    __host__ __device__ auto reciprocal(const T& x) requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        constexpr int GradOrder = T::GradType::Order;
        const auto v = reciprocal(x.template mask<GradOrder>());
        return ResultType(v.value(), -x.grad() * square(v));
    }

    template<Scalar T>
    __host__ __device__ auto sqrt(const T& x) requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        using Tv = T::ValueType;
        constexpr int GradOrder = T::GradType::Order;
        assert(!x.isNegative());
        const auto v = sqrt(x.template mask<GradOrder>());
        return ResultType(v.value(), Tv(0.5) * x.grad() / v);
    }

    template<Scalar T>
    __host__ __device__ auto cbrt(const T& x) requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        using Tv = T::ValueType;
        constexpr int GradOrder = T::GradType::Order;
        const auto v = cbrt(x.template mask<GradOrder>());
        return ResultType(v.value(), Tv(1.0 / 3) * v * x.grad() / x.template mask<GradOrder>());
    }

    template<Scalar T>
    __host__ __device__ auto ln(const T& x) requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        constexpr int GradOrder = T::GradType::Order;
        if constexpr (!T::isComplex)
            assert(x.isPositive() && "[Error]: Invalid param");
        return ResultType(ln(x.value()), x.grad() / x.template mask<GradOrder>());
    }

    template<Scalar T>
    auto ln1p(const T& x) requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        using Tv = T::ValueType;
        constexpr int GradOrder = T::GradType::Order;
        return ResultType(ln1p(x.value()), x.grad() / (Tv(1) + x.template mask<GradOrder>()));
    }

    template<Scalar T, int Order>
    auto ln1pexp(const Diff<T, DiffMode::Forward, Order>& x) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        using GradType = ResultType::GradType;
        constexpr int GradOrder = GradType::Order;
        const GradType x1 = x.template mask<GradOrder>();
        GradType g;
        if (x.value().real().isPositive())
            g = reciprocal(GradType(1) + exp(-x1));
        else {
            const GradType expx = exp(x1);
            g = expx / (GradType(1) + expx);
        }
        return ResultType(ln1pexp(x.value()), g * x.grad());
    }

    //template<Scalar T, int Order>
    //auto log(const Diff<T, DiffMode::Forward, Order>& x, const Diff<T, DiffMode::Forward, Order>& a);

    template<Scalar T>
    __host__ __device__ auto exp(const T& x) requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        constexpr int GradOrder = T::GradType::Order;
        const auto v = exp(x.template mask<GradOrder>());
        return ResultType(v.value(), v * x.grad());
    }

    template<Scalar T, Scalar U>
    __host__ __device__ auto pow(T&& x, U&& a) requires(ForwardDiff<T> && !Diffable<U>) {
        using ResultType = T::ScalarType;
        constexpr int GradOrder = T::GradType::Order;
        const auto y = pow(x.template mask<GradOrder>(), a);
        return ResultType(y.value(), x.grad() * y / x.template mask<GradOrder>() * a);
    }

    //template<Scalar T, int Order>
    //auto pow(const Diff<T, DiffMode::Forward, Order>& x, const Diff<T, DiffMode::Forward, Order>& n);

    template<Scalar T>
    __host__ __device__ auto cos(const T& x) requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        using GradType = T::GradType;
        constexpr int GradOrder = GradType::Order;
        GradType sin_value, cos_value;
        sincos(x.template mask<GradOrder>(), sin_value, cos_value);
        return ResultType(cos_value.value(), -sin_value * x.grad());
    }

    template<Scalar T>
    __host__ __device__ auto sin(const T& x) requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        using GradType = T::GradType;
        constexpr int GradOrder = GradType::Order;
        GradType sin_value, cos_value;
        sincos(x.template mask<GradOrder>(), sin_value, cos_value);
        return ResultType(sin_value.value(), cos_value * x.grad());
    }

    template<Scalar T, Scalar U>
    __host__ __device__ void sincos(const T& x, U& sin_result, U& cos_result) requires(ForwardDiff<T> && ForwardDiff<U>) {
        using ResultType = U::ScalarType;
        using GradType = T::GradType;
        constexpr int GradOrder = GradType::Order;
        GradType sin_value, cos_value;
        sincos(x.template mask<GradOrder>(), sin_value, cos_value);
        sin_result = ResultType(sin_value, cos_value * x.grad());
        cos_result = ResultType(cos_value, -sin_value * x.grad());
    }

    template<Scalar T, int Order>
    auto tan(const Diff<T, DiffMode::Forward, Order>& x) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        constexpr int GradOrder = T::GradType::Order;
        return ResultType(tan(x.value()), x.grad() * square(sec(x.template mask<GradOrder>())));
    }

    template<Scalar T, int Order>
    auto sec(const Diff<T, DiffMode::Forward, Order>& x) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        constexpr int GradOrder = T::GradType::Order;
        const auto x1 = x.template mask<GradOrder>();
        const auto v = sec(x1);
        return ResultType(v.value(), x.grad() * v * tan(x1));
    }

    template<Scalar T, int Order>
    auto csc(const Diff<T, DiffMode::Forward, Order>& x) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        constexpr int GradOrder = T::GradType::Order;
        const auto x1 = x.template mask<GradOrder>();
        const auto v = csc(x1);
        return ResultType(v.value(), -x.grad() * v * cot(x1));
    }

    template<Scalar T, int Order>
    auto cot(const Diff<T, DiffMode::Forward, Order>& x) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        constexpr int GradOrder = T::GradType::Order;
        return ResultType(cot(x.value()), -x.grad() * square(csc(x.template mask<GradOrder>())));
    }

    template<Scalar T, int Order>
    auto arccos(const Diff<T, DiffMode::Forward, Order>& x) {
        using ResultType = Diff<T, DiffMode::Forward, Order>;
        constexpr int GradOrder = T::GradType::Order;
        return ResultType(arccos(x.value()), -x.grad() / sqrt(T(1) - square(x.template mask<GradOrder>())));
    }

    /*
    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<Scalar T, int Order>
    auto arcsin(const Diff<T, DiffMode::Forward, Order>& x);
    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<Scalar T, int Order>
    auto arctan(const Diff<T, DiffMode::Forward, Order>& x);

    template<Scalar T, int Order>
    auto arcsec(const Diff<T, DiffMode::Forward, Order>& x);

    template<Scalar T, int Order>
    auto arccsc(const Diff<T, DiffMode::Forward, Order>& x);

    template<Scalar T, int Order>
    auto arccot(const Diff<T, DiffMode::Forward, Order>& x);
    */

    template<Scalar T>
    __host__ __device__ auto cosh(const T& x) requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        constexpr int GradOrder = T::GradType::Order;
        return ResultType(cosh(x.value()), x.grad() * sinh(x.template mask<GradOrder>()));
    }

    template<Scalar T>
    __host__ __device__ auto sinh(const T& x) requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        constexpr int GradOrder = T::GradType::Order;
        return ResultType(sinh(x.value()), x.grad() * cosh(x.template mask<GradOrder>()));
    }

    template<Scalar T>
    __host__ __device__ auto tanh(const T& x) requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        using Tv = T::ValueType;
        constexpr int GradOrder = T::GradType::Order;
        const auto v = tanh(x.template mask<GradOrder>());
        return ResultType(v.value(), (Tv(1) - square(v)) * x.grad());
    }

    /*
    template<Scalar T, int Order>
    auto sech(const Diff<T, DiffMode::Forward, Order>& x);

    template<Scalar T, int Order>
    auto csch(const Diff<T, DiffMode::Forward, Order>& x);

    template<Scalar T, int Order>
    auto coth(const Diff<T, DiffMode::Forward, Order>& x);

    template<Scalar T, int Order>
    auto arccosh(const Diff<T, DiffMode::Forward, Order>& x);

    template<Scalar T, int Order>
    auto arcsinh(const Diff<T, DiffMode::Forward, Order>& x);

    template<Scalar T, int Order>
    auto arctanh(const Diff<T, DiffMode::Forward, Order>& x);

    template<Scalar T, int Order>
    auto arcsech(const Diff<T, DiffMode::Forward, Order>& x);

    template<Scalar T, int Order>
    auto arccsch(const Diff<T, DiffMode::Forward, Order>& x);

    template<Scalar T, int Order>
    auto arccoth(const Diff<T, DiffMode::Forward, Order>& x);
    */

    template<Scalar T>
    __host__ __device__ auto lncosh(const T& x) requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        constexpr int GradOrder = T::GradType::Order;
        return ResultType(lncosh(x.value()), x.grad() * tanh(x.template mask<GradOrder>()));
    }

    template<Scalar T, int Order>
    T floor(const Diff<T, DiffMode::Forward, Order>& x) {
        return floor(x.value());
    }
    
    template<Scalar T, int Order>
    T ceil(const Diff<T, DiffMode::Forward, Order>& x) {
        return ceil(x.value());
    }
}
