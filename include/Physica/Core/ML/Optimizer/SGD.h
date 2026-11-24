/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "OptBase.h"

namespace Physica {
    /**
     * Stochastic gradient descent for auto diff
     */
    template<Scalar T>
    class SGD : public OptBase<T> {
        T lr;
    public:
        SGD(T lr_);
        SGD(const SGD&) = default;
        SGD(SGD&&) noexcept = default;
        ~SGD() = default;
        /* Operators */
        SGD& operator=(SGD obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void step(Diffable auto& target) const;
        void step(Diffable auto& target, Diffable auto&... targets) const;
        void clear() const noexcept {}
        void swap(SGD& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] auto&& getLearnRate(this auto&& self) noexcept { return self.lr; }
    };

    template<Scalar T>
    SGD<T>::SGD(T lr_) : lr(lr_) {}

    template<Scalar T>
    void SGD<T>::step(Diffable auto& target) const {
        using Target = decltype(target);
        if constexpr (Scalar<Target>)
            target.value() -= lr * target.grad().value();
        else if constexpr (Vector<Target> || Matrix<Target>)
            target.values() -= lr * target.grads().values();
        else
            target.step(*this);
    }

    template<Scalar T>
    void SGD<T>::step(Diffable auto& target, Diffable auto&... targets) const {
        step(target);
        step(targets...);
    }

    template<Scalar T>
    void SGD<T>::swap(SGD& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        lr.swap(obj.lr);
    }
}
