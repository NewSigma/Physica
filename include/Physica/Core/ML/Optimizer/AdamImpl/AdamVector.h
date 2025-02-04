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
#include "AdamImpl.h"

namespace Physica::Core {
    template<Vector V>
    class AdamImpl<V> : public AdamBase<typename V::ScalarType::ValueType> {
        using T = V::ScalarType;
        using Tv = T::ValueType;
        using This = AdamImpl<V>;
        using Base = AdamBase<Tv>;
        using VectorType = DenseVector<Tv, V::SizeAtCompile>;
        using typename Base::Args;

        Tv beta1t;
        Tv beta2t;
        VectorType m;
        VectorType v;
    public:
        AdamImpl(const Args& args, size_t length);
        AdamImpl(const This&) = default;
        AdamImpl(This&&) noexcept = default;
        ~AdamImpl() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void step(const Args& args, void* pTarget) override final;
        void reset(const Args& args) override final;

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return m.getLength(); }
    };

    template<Vector V>
    AdamImpl<V>::AdamImpl(const Args& args, size_t length) : m(length), v(length) {
        reset(args);
    }

    template<Vector V>
    void AdamImpl<V>::step(const Args& args, void* pTarget) {
        V& target = *reinterpret_cast<V*>(pTarget);
        if (!args.decay.isZero())
            target.grads() += args.decay * target.values();

        const Tv beta1 = args.beta1;
        const Tv beta2 = args.beta2;
        const Tv beta1_ = Tv(1) - beta1;
        const Tv beta2_ = Tv(1) - beta2;
        m = beta1 * m + beta1_ * target.grads();
        v = beta2 * v + beta2_ * target.grads().squaredNorms();
        const Tv alpha = args.lr / (Tv(1) - beta1t) * sqrt(Tv(1) - beta2t);
        target.values() -= alpha * divide(m, sqrt(v) + args.epsilon);
        beta1t *= beta1;
        beta2t *= beta2;
    }

    template<Vector V>
    void AdamImpl<V>::reset(const Args& args) {
        beta1t = args.beta1;
        beta2t = args.beta2;
        m = Tv(0);
        v = Tv(0);
    }

    template<Vector V>
    void AdamImpl<V>::swap(This& __restrict obj) noexcept {
        beta1t.swap(obj.beta1t);
        beta2t.swap(obj.beta2t);
        m.swap(obj.m);
        v.swap(obj.v);
    }
}
