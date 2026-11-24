/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "LayerMixin.h"

namespace Physica {
    template<Scalar T>
    class Dropout : public LayerMixin<Dropout<T>> {
        using This = Dropout<T>;
        using Base = LayerMixin<This>;
        using Tv = T::ValueType;

        Tv p;
    public:
        Dropout(Tv p_ = 0.5);
        Dropout(const This&) = default;
        Dropout(This&&) noexcept = default;
        ~Dropout() = default;
        /* Operator */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<RNG R>
        [[nodiscard]] CoDiff<VectorND<T>> forward(const VectorND<T>& x) const;
        void reverse([[maybe_unused]] const This& __restrict obj) const noexcept {}

        void step([[maybe_unused]] auto& optimizer) {}
        void zero_grad() {}

        void swap(This& __restrict obj) noexcept;
        /* Setters */
        void setProb(Tv p_) noexcept { p = p_; }
    };

    template<Scalar T>
    Dropout<T>::Dropout(Tv p_) : p(p_) {
        assert(!p.isNegative() && (p < Tv(1)));
    }

    template<Scalar T>
    template<RNG R>
    CoDiff<VectorND<T>> Dropout<T>::forward(const VectorND<T>& x) const {
        auto mask = VectorND<Tv>::template random_uniform<R>(x.getLength());
        for (auto& elem : mask)
            elem = (elem > p) ? Tv(1) : Tv(0);

        if constexpr (ReverseDiff<T>) {
            auto& y = co_yield hadamard(x.values(), mask);
            x.reverse(hadamard(y.grads(), mask));
        }
        else
            co_return hadamard(x, mask);
    }

    template<Scalar T>
    void Dropout<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        p.swap(obj.p);
    }
}

namespace Physica {
    template<Scalar T>
    class Traits<Dropout<T>> {
    public:
        using ScalarType = T;
    };
}
