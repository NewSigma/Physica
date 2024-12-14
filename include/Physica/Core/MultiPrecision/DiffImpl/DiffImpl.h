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

#include "../Diff.h"

namespace Physica::Core {
    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order>::Diff(T value_) : value(std::move(value_)), grad(0) {}

    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order>::Diff(T value_, GradType grad_) : value(std::move(value_)), grad(std::move(grad_)) {}

    template<Scalar T, DiffMode Mode, int Order>
    template<Scalar U>
    Diff<T, Mode, Order>::Diff(const U& x) {
        static_assert(T::isComplex || !U::isComplex, "[Error]: Cannot convert a complex to a real");
        if constexpr (Diffable<U>) {
            value = x.getValue();
            grad = x.getGrad().template mask<Order - 1>();
        }
        else {
            value = x;
            zero_grad();
        }
    }

    template<Scalar T, DiffMode Mode, int Order>
    inline bool Diff<T, Mode, Order>::operator==(const This& other) const {
        return value == other.value && grad == other.grad;
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<int MaskOrder>
    auto& Diff<T, Mode, Order>::mask() noexcept {
        using MaskedType = std::conditional<MaskOrder == 0, const T&, const Diff<typename Base::ValueType, Mode, MaskOrder>&>::type;
        using ResultType = std::conditional<std::less<int>{}(MaskOrder, Order), MaskedType, const This&>::type;
        return reinterpret_cast<ResultType>(*this);
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<int MaskOrder>
    const auto& Diff<T, Mode, Order>::mask() const noexcept {
        return const_cast<This&>(*this).mask<MaskOrder>();
    }

    template<Scalar T, DiffMode Mode, int Order>
    T Diff<T, Mode, Order>::reverse(GradType grad_) const noexcept {
        static_assert(Mode == DiffMode::Reverse, "[Error]: Call reverse() of a forward diff scalar is not well defined");
        auto& g = const_cast<GradType&>(getGrad());
        g.getValue() += grad_.getValue();
        if constexpr (Order != 1)
            g.reverse(grad_.getGrad());
        return getValue();
    }

    template<Scalar T, DiffMode Mode, int Order>
    inline void Diff<T, Mode, Order>::zero_grad() {
        grad = 0;
    }

    template<Scalar T, DiffMode Mode, int Order>
    auto Diff<T, Mode, Order>::conjugate() const {
        return This(getValue().conjugate(), getGrad().conjugate());
    }

    template<Scalar T, DiffMode Mode, int Order>
    void Diff<T, Mode, Order>::swap(Diff& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        value.swap(obj.value);
        grad.swap(obj.grad);
    }

    template<Scalar T, DiffMode Mode, int Order>
    void Diff<T, Mode, Order>::swap(ScalarRef<This>&& ref) noexcept {
        assert(ScalarPtr<This>(*this) != ScalarPtr<This>(ref) && "[Error]: Self swap is likely a bug");
        value.swap(ref.getValue());
        grad.swap(ref.getGrad());
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<int GradOrder>
    Diff<T, Mode, Order>::Base::template GradRtnTy<GradOrder>&
    Diff<T, Mode, Order>::getGrad() noexcept {
        static_assert(Order >= GradOrder, "[Error]: Order is not enough to calculate the required grad");
        static_assert(GradOrder > 0, "[Error]: 0 or minus order is not well defined");
        if constexpr (GradOrder == 1)
            return grad;
        else
            return getGrad().template getGrad<GradOrder - 1>();
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<int GradOrder>
    inline const Diff<T, Mode, Order>::Base::template GradRtnTy<GradOrder>&
    Diff<T, Mode, Order>::getGrad() const noexcept {
        return const_cast<This&>(*this).template getGrad<GradOrder>();
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ inline bool Diff<T, Mode, Order>::isFinite() const noexcept {
        return value.isFinite() && grad.isFinite();
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<RandomGenerator R>
    inline auto Diff<T, Mode, Order>::random_uniform() {
        return Diff(T::template random_uniform<R>());
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<RandomGenerator R>
    inline auto Diff<T, Mode, Order>::random_normal() {
        return Diff(T::template random_normal<R>());
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<class Distribution, RandomGenerator R>
    inline auto Diff<T, Mode, Order>::random_any(Distribution& dist) {
        return Diff(T::template random_any<Distribution, R>(dist));
    }

#ifdef PHYSICA_HDF5
    template<Scalar T, DiffMode Mode, int Order>
    const H5::DataType& Diff<T, Mode, Order>::getH5DataType() {
        static const auto instance = std::unique_ptr<H5::DataType>([]() -> H5::DataType* {
            auto* result = new H5::DataType(H5T_COMPOUND, sizeof(This));
            const auto id = result->getId();
            H5Tinsert(id, "Value", HOFFSET(This, value), T::getH5DataType().getId());
            H5Tinsert(id, "Grad", HOFFSET(This, grad), GradType::getH5DataType().getId());
            return result;
        }());
        return *instance;
    }
#endif

    template<Scalar T, Scalar U>
    [[nodiscard]] inline CoDiff<typename Internal::BinaryScalarOpRtnTy<std::remove_cvref_t<T>, std::remove_cvref_t<U>>::Type>
    operator+(T&& x, U&& y) requires(Diffable<T>) {
        if constexpr (ForwardDiff<T>) {
            if constexpr (Diffable<U>)
                co_return {x.getValue() + y.getValue(), x.getGrad() + y.getGrad()};
            else
                co_return {x.getValue() + y.getValue(), x.getGrad()};
        }
        else {
            LazyReverse<T&&> x_ = std::forward<T>(x);
            LazyReverse<U&&> y_ = std::forward<U>(y);
            auto result = co_yield x_.getValue() + y_.getValue();
            auto& g = result.getGrad();
            if (!g.isZero()) {
                x_.reverse(g);
                if constexpr (Diffable<U>)
                    y_.reverse(g);
            }
        }
    }

    template<Scalar T, Scalar U>
    [[nodiscard]] inline CoDiff<typename Internal::BinaryScalarOpRtnTy<std::remove_cvref_t<T>, std::remove_cvref_t<U>>::Type>
    operator-(T&& x, U&& y) requires(Diffable<T>) {
        if constexpr (ForwardDiff<T>) {
            if constexpr (Diffable<U>)
                co_return {x.getValue() - y.getValue(), x.getGrad() - y.getGrad()};
            else
                co_return {x.getValue() - y.getValue(), x.getGrad()};
        }
        else {
            LazyReverse<T&&> x_ = std::forward<T>(x);
            LazyReverse<U&&> y_ = std::forward<U>(y);
            auto result = co_yield x_.getValue() - y_.getValue();
            auto& g = result.getGrad();
            if (!g.isZero()) {
                x_.reverse(g);
                if constexpr (Diffable<U>)
                    y_.reverse(-g);
            }
        }
    }

    template<Scalar T, Scalar U>
    [[nodiscard]] inline CoDiff<typename Internal::BinaryScalarOpRtnTy<std::remove_cvref_t<T>, std::remove_cvref_t<U>>::Type>
    operator*(T&& x, U&& y) requires(Diffable<T>) {
        if constexpr (ForwardDiff<T>) {
            using ResultType = typename Internal::BinaryScalarOpRtnTy<std::remove_cvref_t<T>, std::remove_cvref_t<U>>::Type;
            if constexpr (Diffable<U>) {
                using GradType = ResultType::GradType;
                constexpr int GradOrder = GradType::Order;
                co_return ResultType(x.getValue() * y.getValue(), GradType(y.template mask<GradOrder>() * x.getGrad() + x.template mask<GradOrder>() * y.getGrad()));
            }
            else
                co_return ResultType(x.getValue() * y.getValue(), x.getGrad() * y.getValue());
        }
        else {
            LazyReverse<T&&> x_ = std::forward<T>(x);
            LazyReverse<U&&> y_ = std::forward<U>(y);
            auto result = co_yield x_.getValue() * y_.getValue();
            auto& g = result.getGrad();
            if (!g.isZero()) {
                x_.reverse(y_.getValue() * g);
                if constexpr (Diffable<U>)
                    y_.reverse(x_.getValue() * g);
            }
        }
    }

    template<Scalar T, Scalar U>
    [[nodiscard]] inline CoDiff<typename Internal::BinaryScalarOpRtnTy<std::remove_cvref_t<T>, std::remove_cvref_t<U>>::Type>
    operator/(T&& x, U&& y) requires(Diffable<T>) {
        if constexpr (ForwardDiff<T>) {
            using ResultType = typename Internal::BinaryScalarOpRtnTy<std::remove_cvref_t<T>, std::remove_cvref_t<U>>::Type;
            if constexpr (Diffable<U>) {
                using GradType = ResultType::GradType;
                constexpr int GradOrder = GradType::Order;
                const auto v = reciprocal(y.template mask<GradOrder>());
                co_return ResultType(x.getValue() * v.getValue(), GradType((y.template mask<GradOrder>() * x.getGrad() - x.template mask<GradOrder>() * y.getGrad()) * square(v)));
            }
            else {
                const auto v = reciprocal(y.getValue());
                co_return ResultType(x.getValue() * v, x.getGrad() * v);
            }
        }
        else {
            LazyReverse<T&&> x_ = std::forward<T>(x);
            LazyReverse<U&&> y_ = std::forward<U>(y);
            auto result = co_yield x_.getValue() / y_.getValue();
            auto& g = result.getGrad();
            if (!g.isZero()) {
                const auto factor = reciprocal(y_) * g;
                x_.reverse(factor);
                if constexpr (Diffable<U>)
                    y_.reverse(-factor * x_.getValue());
            }
        }
    }

    template<Scalar T, Scalar U>
    [[nodiscard]] inline auto operator+(const U& x, const T& y) requires(Diffable<T> && !Diffable<U>) {
        return y + x;
    }

    template<Scalar T, Scalar U>
    [[nodiscard]] inline auto operator-(const U& x, const T& y) requires(Diffable<T> && !Diffable<U>) {
        return -(y - x);
    }

    template<Scalar T, Scalar U>
    [[nodiscard]] inline auto operator*(const U& x, const T& y) requires(Diffable<T> && !Diffable<U>) {
        return y * x;
    }

    template<Scalar T, Scalar U>
    [[nodiscard]] inline auto operator/(const U& x, const T& y) requires(Diffable<T> && !Diffable<U>) {
        using ResultType = Internal::BinaryScalarOpRtnTy<T, U>::Type;
        const auto rep = reciprocal(y.getValue());
        return ResultType(y.getValue() * rep, -y.getValue() * y.getGrad() * square(rep));
    }

    template<Scalar T>
    inline CoDiff<T> operator-(T&& x) requires(Diffable<T>) {
        if constexpr (ForwardDiff<T>)
            co_return {-x.getValue(), -x.getGrad()};
        else {
            LazyReverse<T&&> x_ = std::forward<T>(x);
            auto result = co_yield x_.getValue();
            auto& g = result.getGrad();
            if (!g.isZero())
                x_.reverse(-g);
        }
    }
}
