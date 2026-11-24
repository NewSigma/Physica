/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Scalar/Diff.h"

namespace Physica {
    template<Packet T>
    [[nodiscard]] auto unit(const T& x) noexcept requires(ForwardDiff<T>) {
        static_assert(!T::isComplex(), "[Error]: Not implemented");
        using ResultType = T::ValueType;
        return ResultType(unit(x.value()));
    }

    template<Packet T1, Packet T2, Packet T3>
    [[nodiscard]] auto fma(const T1 x, const T2 y, const T3 z) noexcept
            requires((ForwardDiff<T1> || ForwardDiff<T2> || ForwardDiff<T3>)
                 && !(ReverseDiff<T1> || ReverseDiff<T2> || ReverseDiff<T3>)) {
        if constexpr (std::same_as<T1, T2> && std::same_as<T2, T3>) {
            using Grad = T1::GradType;
            auto value = fma(x.value(), y.value(), z.value());
            auto grad1 = fma(Grad(x), y.grad(), z.grad());
            auto grad2 = fma(Grad(y), x.grad(), grad1);
            return T1(std::move(value), std::move(grad2));
        }
        else
            return x * y + z;
    }

    template<Packet T>
    [[nodiscard]] T abs(const T& x) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        return T(abs(x.value()), Grad::select(x.value().isPositive(), x.grad(), -x.grad()));
    }

    template<Packet T>
    [[nodiscard]] T square(const T& x) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        return T(square(x.value()), Grad(x) * x.grad() * Grad(2));
    }

    template<Packet T>
    [[nodiscard]] T reciprocal(const T& x) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        const auto y = reciprocal(Grad(x));
        return T(y.value(), -x.grad() * square(y));
    }

    template<Packet T>
    [[nodiscard]] T relu(const T& x) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        return T(relu(x.value()), Grad::select(x.value().isPositive(), x.grad(), Grad(0)));
    }

    template<Packet T>
    [[nodiscard]] T sqrt(const T& x) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        return T(sqrt(x.value()), Grad(0.5) * x.grad() * reciprocal(sqrt(x.value())));
    }

    template<Packet T>
    [[nodiscard]] T cbrt(const T& x) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        return T(cbrt(x.value()), Grad(1.0 / 3) * x.grad() * reciprocal(square(cbrt(x.value()))));
    }

    template<Packet T>
    [[nodiscard]] T ln(const T& x) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        return T(ln(x.value()), x.grad() / Grad(x));
    }

    template<Packet T>
    [[nodiscard]] T ln1p(const T& x) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        return T(ln1p(x.value()), x.grad() / (Grad(1) + Grad(x)));
    }

    template<Packet T>
    [[nodiscard]] T exp(const T& x) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        const auto y = exp(Grad(x));
        return T(y.value(), y * x.grad());
    }

    template<Packet T>
    [[nodiscard]] T expm1(const T& x) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        return T(expm1(x.value()), exp(Grad(x)) * x.grad());
    }

    template<Packet T, Scalar U>
    [[nodiscard]] T pow(const T& x, const U& n) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        return T(pow(x.value(), n), x.grad() * Grad(n) * pow(x.value(), Grad(n) - 1));
    }

    template<Packet T1, Packet T2>
    [[nodiscard]] T1 pow(const T1& x, const T2& y) noexcept requires(ForwardDiff<T1> && ForwardDiff<T2>) {
        const auto v = pow(x.value(), y.value());
        return T1(v, y.value() * reciprocal(x.value()) * x.grad() + v * ln(x.value()) * y.grad());
    }

    template<Packet T>
    [[nodiscard]] T sin(const T& x) noexcept requires(ForwardDiff<T>) {
        return T(sin(x.value()), cos(x.value()) * x.grad());
    }

    template<Packet T>
    [[nodiscard]] T cos(const T& x) noexcept requires(ForwardDiff<T>) {
        return T(cos(x.value()), -sin(x.value()) * x.grad());
    }

    template<Packet T>
    [[nodiscard]] std::pair<T, T> sincos(const T& x) noexcept requires(ForwardDiff<T>) {
        const auto s = sin(x.value());
        const auto c = cos(x.value());
        return {T(s, c * x.grad()), T(c, -s * x.grad())};
    }

    template<Packet T>
    [[nodiscard]] T tan(const T& x) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        const auto v = tan(x.value());
        return T(v, (Grad(1) + square(v)) * x.grad());
    }

    template<Packet T>
    [[nodiscard]] T sec(const T& x) noexcept requires(ForwardDiff<T>) {
        return reciprocal(cos(x));
    }

    template<Packet T>
    [[nodiscard]] T cosh(const T& x) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        return T(cosh(x.value()), x.grad() * sinh(Grad(x)));
    }

    template<Packet T>
    [[nodiscard]] T sinh(const T& x) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        return T(sinh(x.value()), x.grad() * cosh(Grad(x)));
    }

    template<Packet T>
    [[nodiscard]] T tanh(const T& x) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        const Grad y = tanh(Grad(x));
        return T(y.value(), (Grad(1) - square(y)) * x.grad());
    }

    template<Packet T>
    [[nodiscard]] T sech(const T& x) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        const Grad x1 = Grad(x);
        return T(sech(x.value()), -sech(x1) * tanh(x1) * x.grad());
    }

    template<Packet T>
    [[nodiscard]] T arcsinh(const T& x) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        const Grad x1 = Grad(x);
        return T(arcsinh(x.value()), x.grad() / sqrt(square(x1) + Grad(1)));
    }

    template<Packet T>
    [[nodiscard]] T arctanh(const T& x) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        const Grad x1 = Grad(x);
        return T(arctanh(x.value()), x.grad() / (Grad(1) - square(x1)));
    }

    template<Packet T>
    [[nodiscard]] T lncosh(const T& x) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        return T(lncosh(x.value()), tanh(Grad(x)) * x.grad());
    }

    template<Packet T>
    [[nodiscard]] T sigmoid(const T& x) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        const Grad y = sigmoid(Grad(x));
        return T(sigmoid(x.value()), (y - square(y)) * x.grad());
    }

    template<Packet T>
    [[nodiscard]] T softplus(const T& x) noexcept requires(ForwardDiff<T>) {
        return T(softplus(x.value()), sigmoid(x.value()) * x.grad());
    }
}
