/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "../DiffDenseTensor.h"

namespace Physica {
    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    DenseTensor<Diff<T, Mode, Order>, Dims...>::DenseTensor(IndexType shape) : v(std::move(shape)) {
        if constexpr (isForwardDiff())
            g.resize(v);
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    DenseTensor<Diff<T, Mode, Order>, Dims...>::DenseTensor(IndexType shape, T init)
            : v(std::move(shape), init), g(GradTensorType::zeros(v.getShape())) {}

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    DenseTensor<Diff<T, Mode, Order>, Dims...>::DenseTensor(IndexType shape, ScalarType init) requires(isForwardDiff())
            : v(std::move(shape), init.value()), g(v.getShape(), init.grad()) {}

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    DenseTensor<Diff<T, Mode, Order>, Dims...>::DenseTensor(ValueTensor v_, GradTensorType g_)
            : v(std::move(v_)), g(std::move(g_)) {}

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    DenseTensor<Diff<T, Mode, Order>, Dims...>::DenseTensor(const Tensor auto& x) : DenseTensor(x.getShape()) {
        using X = decltype(x);
        static_assert(!ReverseDiff<X>, "[Error]: Copying a reverse diff object breaks the compute graph");
        if constexpr (!Diffable<X>) {
            x.assign(v);
            zero_grad();
        }
        else
            x.assign(*this);
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    DenseTensor<Diff<T, Mode, Order>, Dims...>::DenseTensor(std::same_as<size_t> auto... dims) requires(sizeof...(dims) > 0)
            : DenseTensor(IndexType({dims...})) {}

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    void DenseTensor<Diff<T, Mode, Order>, Dims...>::resize(this auto& self, IndexType shape) {
        self.v.resize(std::move(shape));
        if constexpr (isForwardDiff())
            self.g.resize(self.v);
        else
            assert(self.g.empty() && "[Error]: Reject resize-after-reverse");
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    void DenseTensor<Diff<T, Mode, Order>, Dims...>::reserve(this auto& self, size_t size) noexcept {
        self.v.reserve(size);
        if constexpr (isForwardDiff())
            self.g.reserve(size);
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    template<RNG R>
    void DenseTensor<Diff<T, Mode, Order>, Dims...>::random_uniform(this auto& self) {
        self.v.template random_uniform<R>();
        self.zero_grad();
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    template<RNG R>
    void DenseTensor<Diff<T, Mode, Order>, Dims...>::random_normal(this auto& self) {
        self.v.template random_normal<R>();
        self.zero_grad();
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    template<RNG R>
    void DenseTensor<Diff<T, Mode, Order>, Dims...>::random_any(this auto& self, auto& distribution) {
        self.v.template random_any<R>(distribution);
        self.zero_grad();
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    void DenseTensor<Diff<T, Mode, Order>, Dims...>::zero_grad(this auto& self) {
        self.g.zeros();
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    void DenseTensor<Diff<T, Mode, Order>, Dims...>::junk(this auto& self) {
        self.v.junk();
        if constexpr (isForwardDiff())
            self.g.junk();
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    void DenseTensor<Diff<T, Mode, Order>, Dims...>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        v.swap(obj.v);
        g.swap(obj.g);
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    auto DenseTensor<Diff<T, Mode, Order>, Dims...>::data_handle(this auto&& self) noexcept {
        constexpr bool IsConst = std::is_const<std::remove_reference_t<decltype(self)>>::value;
        using U = Diff<T, Mode, Order>;
        using RetTy = std::conditional<IsConst, typename U::ConstPtrTy, typename U::PtrTy>::type;
        if constexpr (isReverseDiff())
            assert(!self.g.empty() && "[Error]: Grad is not ready");
        return RetTy(self.v.data(), self.g.data());
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    auto&& DenseTensor<Diff<T, Mode, Order>, Dims...>::values(this auto&& self) noexcept {
        return forward_like<decltype(self)>(self.v);
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    template<int GradOrder>
    auto&& DenseTensor<Diff<T, Mode, Order>, Dims...>::grads(this auto&& self) noexcept {
        if constexpr (isReverseDiff()) {
            if (self.g.empty()) {
                auto& g1 = const_cast<GradTensorType&>(self.g);
                g1.resize(self.v);
                g1.zeros();
            }
        }

        if constexpr (GradOrder == 1)
            return forward_like<decltype(self)>(self.g);
        else
            return forward_like<decltype(self)>(self.g.template grads<GradOrder - 1>());
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    auto DenseTensor<Diff<T, Mode, Order>, Dims...>::zeros(IndexType shape) -> This {
        return This(ValueTensor::zeros(std::move(shape)));
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    auto DenseTensor<Diff<T, Mode, Order>, Dims...>::zeros(std::same_as<size_t> auto... dims) -> This requires(sizeof...(dims) > 0) {
        static_assert(sizeof...(dims) == This::ndim(), "[Error]: NDim is not consistent");
        return zeros(IndexType({dims...}));
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    auto DenseTensor<Diff<T, Mode, Order>, Dims...>::junk(IndexType shape) -> This {
        This result(std::move(shape));
        result.junk();
        return result;
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    auto DenseTensor<Diff<T, Mode, Order>, Dims...>::junk(std::same_as<size_t> auto... dims) -> This requires(sizeof...(dims) > 0) {
        static_assert(sizeof...(dims) == This::ndim(), "[Error]: NDim is not consistent");
        return junk(IndexType({dims...}));
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    template<RNG R>
    auto DenseTensor<Diff<T, Mode, Order>, Dims...>::random_uniform(IndexType shape) -> This {
        return This(ValueTensor::template random_uniform<R>(std::move(shape)));
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    template<RNG R>
    auto DenseTensor<Diff<T, Mode, Order>, Dims...>::random_uniform(std::same_as<size_t> auto... dims) -> This requires(sizeof...(dims) > 0) {
        static_assert(sizeof...(dims) == This::ndim(), "[Error]: NDim is not consistent");
        return random_uniform<R>(IndexType({dims...}));
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    template<RNG R>
    auto DenseTensor<Diff<T, Mode, Order>, Dims...>::random_normal(IndexType shape) -> This {
        return This(ValueTensor::template random_normal<R>(std::move(shape)));
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    template<RNG R>
    auto DenseTensor<Diff<T, Mode, Order>, Dims...>::random_normal(std::same_as<size_t> auto... dims) -> This requires(sizeof...(dims) > 0) {
        static_assert(sizeof...(dims) == This::ndim(), "[Error]: NDim is not consistent");
        return random_normal<R>(IndexType({dims...}));
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    template<RNG R>
    auto DenseTensor<Diff<T, Mode, Order>, Dims...>::random_any(IndexType shape, auto& distribution) -> This {
        return This(ValueTensor::template random_any<R>(std::move(shape), distribution));
    }
}
