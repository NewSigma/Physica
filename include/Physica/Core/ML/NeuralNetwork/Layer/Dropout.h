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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "LayerBase.h"

namespace Physica {
    template<Scalar T>
    class Dropout : public LayerBase<Dropout<T>> {
        using This = Dropout<T>;
        using Base = LayerBase<This>;
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
        void reverse(const This& __restrict obj) const noexcept {}

        void step(auto& optimizer) {}
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
            auto y = co_yield hadamard(x.values(), mask);
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
