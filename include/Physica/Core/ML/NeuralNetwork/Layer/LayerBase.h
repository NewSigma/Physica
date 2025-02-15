/*
 * Copyright 2023-2025 Weibo He.
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

#include <concepts>
#include "Physica/CRTPBase.h"

namespace Physica {
    template<class Derived>
    class LayerBase : public CRTPBase<LayerBase<Derived>> {
        using This = LayerBase<Derived>;
        using Base = CRTPBase<This>;
        using TraitsType = Traits<Derived>;
    public:
        using ScalarType = TraitsType::ScalarType;
        constexpr static bool IsTrain = ScalarType::isDiffable;
        constexpr static bool IsInfer = !IsTrain;
    public:
        ~LayerBase() = default;
        /* Operations */
        template<class T>
        [[nodiscard]] auto forward(const T& x) const { return Base::getDerived().template forward<T>(x); }
        auto reverse(const Derived& __restrict other) const noexcept;

        template<class Optimizer>
        auto step(Optimizer& opt) { return Base::getDerived().step(opt); }
        auto step() { return Base::getDerived().step(); }
        auto zero_grad() { return Base::getDerived().zero_grad(); }
    protected:
        LayerBase() = default;
        LayerBase(const This&) = default;
        LayerBase(This&&) noexcept = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
    };

    template<class Derived>
    auto LayerBase<Derived>::reverse(const Derived& __restrict other) const noexcept {
        assert(this != &other && "[Error]: Self reverse is invalid");
        return Base::getDerived().reverse(other);
    }
    /**
     * Deep Neutral Network
     */
    template<class T>
    concept DNN = std::derived_from<T, LayerBase<T>>;
}

namespace Physica {
    template<class T>
    class Traits<LayerBase<T>> {
    public:
        using Derived = T;
    };
}
