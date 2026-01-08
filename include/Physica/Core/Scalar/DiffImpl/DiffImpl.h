/*
 * Copyright 2024-2026 Weibo He.
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

#include "../Diff.h"

namespace Physica {
    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ Diff<T, Mode, Order>::Diff(T v_) : v(std::move(v_)), g(0) {}

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ Diff<T, Mode, Order>::Diff(T v_, GradType g_) : v(std::move(v_)), g(std::move(g_)) {}

    template<Scalar T, DiffMode Mode, int Order>
    template<Scalar U>
    __host__ __device__ Diff<T, Mode, Order>::Diff(const U& x) {
        Base::template assert_assign<U>();
        if constexpr (Diffable<U>) {
            v = x.value();
            g = x.grad().template mask<Order - 1>();
        }
        else {
            v = x;
            zero_grad();
        }
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ bool Diff<T, Mode, Order>::operator==(const This& other) const {
        return v == other.v && g == other.g;
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<int MaskOrder>
    __host__ __device__ auto Diff<T, Mode, Order>::mask() const noexcept {
        using MaskedType = std::conditional<MaskOrder == 0, T, Diff<typename Base::ValueType, Mode, MaskOrder>>::type;
        using ResultType = std::conditional<(MaskOrder < Order), MaskedType, This>::type;
        return reinterpret_cast<const ResultType&>(*this);
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ T Diff<T, Mode, Order>::reverse(GradType grad) const noexcept {
        static_assert(Mode == DiffMode::Reverse, "[Error]: Call reverse() of a forward diff scalar is not well defined");
        auto& g1 = const_cast<GradType&>(this->grad());
        g1.value() += grad.value();
        if constexpr (Order != 1)
            g1.reverse(grad.grad());
        return v;
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ void Diff<T, Mode, Order>::zero_grad() {
        g = 0;
    }

    template<Scalar T, DiffMode Mode, int Order>
    auto Diff<T, Mode, Order>::conjugate() const noexcept -> This {
        return This(v.conjugate(), g.conjugate());
    }

    template<Scalar T, DiffMode Mode, int Order>
    auto Diff<T, Mode, Order>::squaredNorm() const noexcept -> CoDiff<This> {
        if constexpr (isReverseDiff) {
            auto y = co_yield value().squaredNorm();
            reverse(T(2) * v * y.grad());
        }
        else
            co_return This(v.squaredNorm(), T(2) * v * g);
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ void Diff<T, Mode, Order>::swap(Diff& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        v.swap(obj.v);
        g.swap(obj.g);
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ auto Diff<T, Mode, Order>::real_ptr() noexcept {
        return RealType::PtrTy(v.real_ptr(), g.real_ptr());
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ const auto Diff<T, Mode, Order>::real_ptr() const noexcept {
        return Base::real_ptr();
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ void Diff<T, Mode, Order>::swap(ScalarRef<This>&& ref) noexcept {
        assert(ScalarPtr<This>(*this) != ScalarPtr<This>(ref) && "[Error]: Self swap is likely a bug");
        v.swap(ref.value());
        g.swap(ref.grad());
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<int GradOrder>
    __host__ __device__ auto& Diff<T, Mode, Order>::grad() noexcept {
        static_assert(Order >= GradOrder, "[Error]: Order is not enough to calculate the required grad");
        static_assert(GradOrder > 0, "[Error]: 0 or minus order is not well defined");
        if constexpr (GradOrder == 1)
            return g;
        else
            return grad().template grad<GradOrder - 1>();
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<int GradOrder>
    __host__ __device__ const auto& Diff<T, Mode, Order>::grad() const noexcept {
        return const_cast<This&>(*this).template grad<GradOrder>();
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ bool Diff<T, Mode, Order>::isFinite() const noexcept {
        return v.isFinite() && g.isFinite();
    }

    template<Scalar T, DiffMode Mode, int Order>
    auto Diff<T, Mode, Order>::fromPhase(RealType phase) noexcept -> This requires(isComplex) {
        static_assert(isForwardDiff && (Order == 1), "[Error]: Not implemented");
        T y = T::fromPhase(phase.value());
        return This(y, T(-y.imag(), y.real()) * phase.grad());
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<RNG R>
    auto Diff<T, Mode, Order>::random_uniform() {
        return Diff(T::template random_uniform<R>());
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<RNG R>
    auto Diff<T, Mode, Order>::random_normal() {
        return Diff(T::template random_normal<R>());
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<RNG R>
    auto Diff<T, Mode, Order>::random_any(auto& distribution) {
        return Diff(T::template random_any<R>(distribution));
    }

#ifdef PHYSICA_HDF5
    template<Scalar T, DiffMode Mode, int Order>
    const H5::DataType& Diff<T, Mode, Order>::dtype_hdf5() noexcept {
        static const auto instance = std::unique_ptr<H5::DataType>([]() -> H5::DataType* {
            auto* result = new (std::nothrow) H5::DataType(H5T_COMPOUND, sizeof(This));
            const auto id = result->getId();
            H5Tinsert(id, "Value", HOFFSET(This, v), T::dtype_hdf5().getId());
            H5Tinsert(id, "Grad", HOFFSET(This, g), GradType::dtype_hdf5().getId());
            return result;
        }());
        return *instance;
    }
#endif

    template<Scalar T, Scalar U>
    [[nodiscard]] __host__ __device__ auto operator+(const T& x, const U& y) noexcept requires(ForwardDiff<T>) {
        using ResultType = Internal::BinaryScalarOpRtnTy<T, U>::Type;
        if constexpr (Diffable<U>)
            return ResultType{x.value() + y.value(), x.grad() + y.grad()};
        else
            return ResultType{x.value() + y.value(), x.grad()};
    }

    template<Scalar T, Scalar U>
    [[nodiscard]] auto operator+(T&& x, U&& y) noexcept -> CoDiff<typename Internal::BinaryScalarOpRtnTy<std::remove_cvref_t<T>, std::remove_cvref_t<U>>::Type>
            requires(ReverseDiff<T>) {
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        LazyDestroy<U&&> y_ = std::forward<U>(y);
        auto result = co_yield x_.value() + y_.value();
        auto& g = result.grad();
        x_.reverse(g);
        if constexpr (Diffable<U>)
            y_.reverse(g);
    }

    template<Scalar T, Scalar U>
    [[nodiscard]] __host__ __device__ auto operator-(const T& x, const U& y) noexcept requires(ForwardDiff<T>) {
        using ResultType = Internal::BinaryScalarOpRtnTy<T, U>::Type;
        if constexpr (Diffable<U>)
            return ResultType{x.value() - y.value(), x.grad() - y.grad()};
        else
            return ResultType{x.value() - y.value(), x.grad()};
    }

    template<Scalar T, Scalar U>
    [[nodiscard]] auto operator-(T&& x, U&& y) noexcept -> CoDiff<typename Internal::BinaryScalarOpRtnTy<std::remove_cvref_t<T>, std::remove_cvref_t<U>>::Type>
            requires(ReverseDiff<T>) {
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        LazyDestroy<U&&> y_ = std::forward<U>(y);
        auto result = co_yield x_.value() - y_.value();
        auto& g = result.grad();
        x_.reverse(g);
        if constexpr (Diffable<U>)
            y_.reverse(-g);
    }

    template<Scalar T, Scalar U>
    [[nodiscard]] __host__ __device__ auto operator*(const T& x, const U& y) noexcept requires(ForwardDiff<T>) {
        using ResultType = typename Internal::BinaryScalarOpRtnTy<T, U>::Type;
        if constexpr (Diffable<U>) {
            using GradType = ResultType::GradType;
            constexpr int GradOrder = GradType::Order;
            return ResultType(x.value() * y.value(), GradType(y.template mask<GradOrder>() * x.grad() + x.template mask<GradOrder>() * y.grad()));
        }
        else
            return ResultType(x.value() * y.value(), x.grad() * y.value());
    }

    template<Scalar T, Scalar U>
    [[nodiscard]] auto operator*(T&& x, U&& y) noexcept -> CoDiff<typename Internal::BinaryScalarOpRtnTy<std::remove_cvref_t<T>, std::remove_cvref_t<U>>::Type>
            requires(ReverseDiff<T>) {
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        LazyDestroy<U&&> y_ = std::forward<U>(y);
        auto result = co_yield x_.value() * y_.value();
        auto& g = result.grad();
        x_.reverse(y_.value() * g);
        if constexpr (Diffable<U>)
            y_.reverse(x_.value() * g);
    }

    template<Scalar T, Scalar U>
    [[nodiscard]] __host__ __device__ auto operator/(const T& x, const U& y) noexcept requires(ForwardDiff<T>) {
        using ResultType = typename Internal::BinaryScalarOpRtnTy<std::remove_cvref_t<T>, std::remove_cvref_t<U>>::Type;
        if constexpr (Diffable<U>) {
            using GradType = ResultType::GradType;
            constexpr int GradOrder = GradType::Order;
            const auto v = reciprocal(y.template mask<GradOrder>());
            return ResultType(x.value() * v.value(), GradType((y.template mask<GradOrder>() * x.grad() - x.template mask<GradOrder>() * y.grad()) * square(v)));
        }
        else {
            const auto v = reciprocal(y.value());
            return ResultType(x.value() * v, x.grad() * v);
        }
    }

    template<Scalar T, Scalar U>
    [[nodiscard]] auto operator/(T&& x, U&& y) noexcept -> CoDiff<typename Internal::BinaryScalarOpRtnTy<std::remove_cvref_t<T>, std::remove_cvref_t<U>>::Type>
            requires(ReverseDiff<T>) {
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        LazyDestroy<U&&> y_ = std::forward<U>(y);
        auto result = co_yield x_.value() / y_.value();
        const auto factor = reciprocal(y_.value()) * result.grad();
        x_.reverse(factor);
        if constexpr (Diffable<U>)
            y_.reverse(-factor * result.value());
    }

    template<Scalar T, Scalar U>
    [[nodiscard]] __host__ __device__ auto operator+(U&& x, T&& y) noexcept requires(Diffable<T> && !Diffable<U>) {
        if constexpr (IsHost() || ForwardDiff<T>)
            return std::forward<T>(y) + std::forward<U>(x);
        else
            unreachable();
    }

    template<Scalar T, Scalar U>
    [[nodiscard]] __host__ __device__ auto operator-(U&& x, T&& y) noexcept requires(Diffable<T> && !Diffable<U>) {
        if constexpr (IsHost() || ForwardDiff<T>)
            return -(std::forward<T>(y) - std::forward<U>(x));
        else
            unreachable();
    }

    template<Scalar T, Scalar U>
    [[nodiscard]] __host__ __device__ auto operator*(U&& x, T&& y) noexcept requires(Diffable<T> && !Diffable<U>) {
        if constexpr (IsHost() || ForwardDiff<T>)
            return std::forward<T>(y) * std::forward<U>(x);
        else
            unreachable();
    }

    template<Scalar T, Scalar U>
    [[nodiscard]] __host__ __device__ auto operator/(const U& x, const T& y) noexcept requires(ForwardDiff<T> && !Diffable<U>) {
        using ResultType = Internal::BinaryScalarOpRtnTy<T, U>::Type;
        const auto rep = reciprocal(y.value());
        return ResultType(x * rep, -x * square(rep) * y.grad());
    }

    template<Scalar T, Scalar U>
    [[nodiscard]] auto operator/(U&& x, T&& y) noexcept -> CoDiff<typename Internal::BinaryScalarOpRtnTy<std::remove_cvref_t<T>, std::remove_cvref_t<U>>::Type>
            requires(ReverseDiff<T> && !Diffable<U>) {
        const auto rep = reciprocal(std::forward<T>(y));
        auto result = co_yield x * rep.value();
        rep.reverse(x * result.grad());
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ auto operator-(const T& x) noexcept requires(ForwardDiff<T>) {
        using ResultType = T::ScalarType;
        return ResultType(-x.value(), -x.grad());
    }

    template<Scalar T>
    [[nodiscard]] CoDiff<T> operator-(T&& x) noexcept requires(ReverseDiff<T>) {
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        auto y = co_yield -x_.value();
        x_.reverse(-y.grad());
    }
}
