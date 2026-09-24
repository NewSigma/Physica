/*
 * Copyright 2025-2026 Weibo He.
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
    [[nodiscard]] __host__ __device__ auto unit(const T& x) noexcept requires(ForwardDiff<T>) {
        static_assert(!T::isComplex(), "[Error]: Not implemented");
        using ResultType = T::ValueType;
        return ResultType(unit(x.value()));
    }

    template<Scalar T1, Scalar T2, Scalar T3>
    [[nodiscard]] __host__ __device__ auto fma(const T1 x, const T2 y, const T3 z) noexcept
            requires((ForwardDiff<T1> || ForwardDiff<T2> || ForwardDiff<T3>)
                 && !(ReverseDiff<T1> || ReverseDiff<T2> || ReverseDiff<T3>)) {
        if constexpr (std::same_as<T1, T2> && std::same_as<T2, T3>) {
            constexpr int Order = T1::Order;
            const auto value = fma(x.value(), y.value(), z.value());
            const auto grad1 = fma(x.template grad_mask<Order - 1>(), y.grad(), z.grad());
            const auto grad2 = fma(y.template grad_mask<Order - 1>(), x.grad(), grad1);
            return T1(std::move(value), std::move(grad2));
        }
        else
            return x * y + z;
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ auto abs(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        return ResultType(abs(x.value()), x.value().isPositive() ? x.grad() : -x.grad());
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ auto relu(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        using GradType = T::GradType;
        return ResultType(relu(x.value()), x.value().isPositive() ? x.grad() : GradType(0));
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ auto square(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        using Tv = T::ValueType;
        constexpr int GradOrder = T::GradType::Order;
        return ResultType(square(x.value()), Tv(2) * x.template grad_mask<GradOrder>() * x.grad());
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ auto reciprocal(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        constexpr int GradOrder = T::GradType::Order;
        const auto v = reciprocal(x.template grad_mask<GradOrder>());
        return ResultType(v.value(), -x.grad() * square(v));
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ auto sqrt(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        using Tv = T::ValueType;
        constexpr int GradOrder = T::GradType::Order;
        if constexpr (!T::isComplex())
            assert(!x.isNegative());
        const auto v = sqrt(x.template grad_mask<GradOrder>());
        return ResultType(v.value(), Tv(0.5) * x.grad() / v);
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ auto cbrt(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        using Tv = T::ValueType;
        constexpr int GradOrder = T::GradType::Order;
        const auto v = cbrt(x.template grad_mask<GradOrder>());
        return ResultType(v.value(), Tv(1.0 / 3) * v * x.grad() / x.template grad_mask<GradOrder>());
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ auto ln(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        constexpr int GradOrder = T::GradType::Order;
        if constexpr (!T::isComplex())
            assert(x.isPositive() && "[Error]: Invalid param");
        return ResultType(ln(x.value()), x.grad() / x.template grad_mask<GradOrder>());
    }

    template<Scalar T>
    [[nodiscard]] auto ln1p(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        using Tv = T::ValueType;
        constexpr int GradOrder = T::GradType::Order;
        return ResultType(ln1p(x.value()), x.grad() / (Tv(1) + x.template grad_mask<GradOrder>()));
    }

    template<Scalar T, Scalar U>
    [[nodiscard]] auto log(const T& x, const U& a) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        constexpr int GradOrder = T::GradType::Order;
        const auto x1 = x.template grad_mask<GradOrder>();
        const auto lna = ln(a.value());
        auto grad = x.grad() / (x1 * lna);
        if constexpr (Diffable<U>)
            grad = grad - ln(x1) * a.grad() / (a.template grad_mask<GradOrder>() * square(lna));
        return ResultType(log(x.value(), a.value()), std::move(grad));
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ auto exp(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        constexpr int GradOrder = T::GradType::Order;
        const auto v = exp(x.template grad_mask<GradOrder>());
        return ResultType(v.value(), v * x.grad());
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ auto expm1(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        constexpr int GradOrder = T::GradType::Order;
        return ResultType(expm1(x.value()), exp(x.template grad_mask<GradOrder>()) * x.grad());
    }

    template<Scalar T, Scalar U>
    [[nodiscard]] __host__ __device__ auto pow(T&& x, U&& a) noexcept requires(ForwardDiff<T> && !Diffable<U>) {
        using ResultType = T::ScalarType;
        constexpr int GradOrder = T::GradType::Order;
        const auto y = pow(x.template grad_mask<GradOrder>(), a);
        return ResultType(y.value(), x.grad() * y / x.template grad_mask<GradOrder>() * a);
    }

    template<Scalar T, Scalar U>
    [[nodiscard]] auto pow(const T& x, const U& n) noexcept requires(ForwardDiff<T> && ForwardDiff<U>) {
        using ResultType = T::ScalarType;
        constexpr int GradOrder = T::GradType::Order;
        const auto x1 = x.template grad_mask<GradOrder>();
        const auto n1 = n.template grad_mask<GradOrder>();
        const auto y = pow(x1, n1);
        return ResultType(y.value(), x.grad() * y / x1 * n1 + y * ln(x1) * n.grad());
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ auto cos(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        using GradType = T::GradType;
        constexpr int GradOrder = GradType::Order;
        auto [s, c] = sincos(x.template grad_mask<GradOrder>());
        return ResultType(c.value(), -s * x.grad());
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ auto sin(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        using GradType = T::GradType;
        constexpr int GradOrder = GradType::Order;
        auto [s, c] = sincos(x.template grad_mask<GradOrder>());
        return ResultType(s.value(), c * x.grad());
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ std::pair<T, T> sincos(const T& x) noexcept requires(ForwardDiff<T>) {
        auto [s, c] = sincos(x.template grad_mask<T::GradType::Order>());
        return {T(s, c * x.grad()), T(c, -s * x.grad())};
    }

    template<Scalar T>
    [[nodiscard]] auto tan(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        constexpr int GradOrder = T::GradType::Order;
        return ResultType(tan(x.value()), x.grad() * square(sec(x.template grad_mask<GradOrder>())));
    }

    template<Scalar T>
    [[nodiscard]] auto sec(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        constexpr int GradOrder = T::GradType::Order;
        const auto x1 = x.template grad_mask<GradOrder>();
        const auto v = sec(x1);
        return ResultType(v.value(), x.grad() * v * tan(x1));
    }

    template<Scalar T>
    [[nodiscard]] auto csc(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        constexpr int GradOrder = T::GradType::Order;
        const auto x1 = x.template grad_mask<GradOrder>();
        const auto v = csc(x1);
        return ResultType(v.value(), -x.grad() * v * cot(x1));
    }

    template<Scalar T>
    [[nodiscard]] auto cot(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        constexpr int GradOrder = T::GradType::Order;
        return ResultType(cot(x.value()), -x.grad() * square(csc(x.template grad_mask<GradOrder>())));
    }

    template<Scalar T>
    [[nodiscard]] auto arccos(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        using Tv = T::ValueType;
        constexpr int GradOrder = T::GradType::Order;
        return ResultType(arccos(x.value()), -x.grad() / sqrt(Tv(1) - square(x.template grad_mask<GradOrder>())));
    }

    template<Scalar T>
    [[nodiscard]] auto arcsin(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        using Tv = T::ValueType;
        constexpr int GradOrder = T::GradType::Order;
        const auto x1 = x.template grad_mask<GradOrder>();
        return ResultType(arcsin(x.value()), x.grad() / sqrt(Tv(1) - square(x1)));
    }


    template<Scalar T>
    [[nodiscard]] auto arctan(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        using Tv = T::ValueType;
        constexpr int GradOrder = T::GradType::Order;
        const auto x1 = x.template grad_mask<GradOrder>();
        return ResultType(arctan(x.value()), x.grad() / (Tv(1) + square(x1)));
    }

    template<Scalar T>
    [[nodiscard]] auto arcsec(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        using Tv = T::ValueType;
        constexpr int GradOrder = T::GradType::Order;
        const auto x1 = x.template grad_mask<GradOrder>();
        const auto u = reciprocal(x1);
        return ResultType(arcsec(x.value()), x.grad() / (square(x1) * sqrt(Tv(1) - square(u))));
    }

    template<Scalar T>
    [[nodiscard]] auto arccsc(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        using Tv = T::ValueType;
        constexpr int GradOrder = T::GradType::Order;
        const auto x1 = x.template grad_mask<GradOrder>();
        const auto u = reciprocal(x1);
        return ResultType(arccsc(x.value()), -x.grad() / (square(x1) * sqrt(Tv(1) - square(u))));
    }

    template<Scalar T>
    [[nodiscard]] auto arccot(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        using Tv = T::ValueType;
        constexpr int GradOrder = T::GradType::Order;
        const auto x1 = x.template grad_mask<GradOrder>();
        return ResultType(arccot(x.value()), -x.grad() / (Tv(1) + square(x1)));
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ auto cosh(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        constexpr int GradOrder = T::GradType::Order;
        return ResultType(cosh(x.value()), x.grad() * sinh(x.template grad_mask<GradOrder>()));
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ auto sinh(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        constexpr int GradOrder = T::GradType::Order;
        return ResultType(sinh(x.value()), x.grad() * cosh(x.template grad_mask<GradOrder>()));
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ auto tanh(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        using Tv = T::ValueType;
        constexpr int GradOrder = T::GradType::Order;
        const auto v = tanh(x.template grad_mask<GradOrder>());
        return ResultType(v.value(), (Tv(1) - square(v)) * x.grad());
    }

    template<Scalar T>
    [[nodiscard]] auto sech(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        constexpr int GradOrder = T::GradType::Order;
        const auto x1 = x.template grad_mask<GradOrder>();
        return ResultType(sech(x.value()), -sech(x1) * tanh(x1) * x.grad());
    }

    template<Scalar T>
    [[nodiscard]] auto csch(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        constexpr int GradOrder = T::GradType::Order;
        const auto x1 = x.template grad_mask<GradOrder>();
        return ResultType(csch(x.value()), -csch(x1) * coth(x1) * x.grad());
    }

    template<Scalar T>
    [[nodiscard]] auto coth(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        constexpr int GradOrder = T::GradType::Order;
        const auto x1 = x.template grad_mask<GradOrder>();
        return ResultType(coth(x.value()), -square(csch(x1)) * x.grad());
    }

    template<Scalar T>
    [[nodiscard]] auto arccosh(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        using Tv = T::ValueType;
        constexpr int GradOrder = T::GradType::Order;
        const auto x1 = x.template grad_mask<GradOrder>();
        return ResultType(arccosh(x.value()), x.grad() / sqrt(square(x1) - Tv(1)));
    }

    template<Scalar T>
    [[nodiscard]] auto arcsinh(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        using Tv = T::ValueType;
        constexpr int GradOrder = T::GradType::Order;
        const auto x1 = x.template grad_mask<GradOrder>();
        return ResultType(arcsinh(x.value()), x.grad() / sqrt(square(x1) + Tv(1)));
    }

    template<Scalar T>
    [[nodiscard]] auto arctanh(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        using Tv = T::ValueType;
        constexpr int GradOrder = T::GradType::Order;
        const auto x1 = x.template grad_mask<GradOrder>();
        return ResultType(arctanh(x.value()), x.grad() / (Tv(1) - square(x1)));
    }

    template<Scalar T>
    [[nodiscard]] auto arcsech(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        using Tv = T::ValueType;
        constexpr int GradOrder = T::GradType::Order;
        const auto x1 = x.template grad_mask<GradOrder>();
        const auto u = reciprocal(x1);
        return ResultType(arcsech(x.value()), -x.grad() / (square(x1) * sqrt(square(u) - Tv(1))));
    }

    template<Scalar T>
    [[nodiscard]] auto arccsch(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        using Tv = T::ValueType;
        constexpr int GradOrder = T::GradType::Order;
        const auto x1 = x.template grad_mask<GradOrder>();
        const auto u = reciprocal(x1);
        return ResultType(arccsch(x.value()), -x.grad() / (square(x1) * sqrt(square(u) + Tv(1))));
    }

    template<Scalar T>
    [[nodiscard]] auto arccoth(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        using Tv = T::ValueType;
        constexpr int GradOrder = T::GradType::Order;
        const auto x1 = x.template grad_mask<GradOrder>();
        return ResultType(arccoth(x.value()), x.grad() / (Tv(1) - square(x1)));
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ auto lncosh(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        constexpr int GradOrder = T::GradType::Order;
        return ResultType(lncosh(x.value()), x.grad() * tanh(x.template grad_mask<GradOrder>()));
    }

    template<Scalar T>
    [[nodiscard]] auto softplus(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        constexpr int GradOrder = T::GradType::Order;
        return ResultType(softplus(x.value()), x.grad() * sigmoid(x.template grad_mask<GradOrder>()));
    }

    template<Scalar T>
    [[nodiscard]] auto floor(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ValueType;
        return ResultType(floor(x.value()));
    }

    template<Scalar T>
    [[nodiscard]] auto ceil(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ValueType;
        return ResultType(ceil(x.value()));
    }
}
