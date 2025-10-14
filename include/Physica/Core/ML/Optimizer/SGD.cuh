/*
 * Copyright 2024-2025 Weibo He.
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

#include "SGD.h"

namespace Physica {
    /**
     * Stochastic gradient descent for auto diff
     */
    template<Scalar T>
    class device_obj<SGD<T>> {
        static_assert(!Diffable<T>);
        using host_obj = SGD<T>;
        using This = device_obj<SGD<T>>;

        T lr;
    public:
        device_obj(T lr_);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void step(Diffable auto& target) const;
        void step(Diffable auto& target, Diffable auto&... targets) const;

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] T& getLearnRate() noexcept { return lr; }
        [[nodiscard]] const T& getLearnRate() const noexcept { return lr; }
    };

    template<Scalar T>
    device_obj<SGD<T>>::device_obj(T lr_) : lr(lr_) {}

    template<Scalar T>
    void device_obj<SGD<T>>::step(Diffable auto& target) const {
        using U = decltype(target);
        if constexpr (Scalar<U>)
            target.value() -= lr * target.grad().value();
        else if constexpr (Vector<U> || Matrix<U>)
            target.values() -= lr * target.grads().values();
        else
            target.step(*this);
    }

    template<Scalar T>
    void device_obj<SGD<T>>::step(Diffable auto& target, Diffable auto&... targets) const {
        step(target);
        step(targets...);
    }

    template<Scalar T>
    void device_obj<SGD<T>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        lr.swap(obj.lr);
    }
}
