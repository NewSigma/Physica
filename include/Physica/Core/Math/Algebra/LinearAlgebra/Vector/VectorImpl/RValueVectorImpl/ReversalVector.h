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

#include "../RValueVector.h"

namespace Physica {
    template<Vector V>
    class ReversalVector final : public RValueVector<ReversalVector<V>> {
        using This = ReversalVector<V>;
        using Base = RValueVector<This>;

        V& v;
    public:
        using typename Base::ScalarType;
    public:
        explicit ReversalVector(V& v_) : v(v_) {}
        ReversalVector(const This&) = delete;
        ReversalVector(This&&) noexcept = delete;
        ~ReversalVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t index) const { return v.calc(getLength() - index - 1); }
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
    };

    template<class Derived, Scalar ScalarT>
    auto RValueVector<Derived, ScalarT>::reversal() noexcept {
        return ReversalVector<Derived>(Base::getDerived());
    }

    template<class Derived, Scalar ScalarT>
    const auto RValueVector<Derived, ScalarT>::reversal() const noexcept {
        return ReversalVector<Derived>(Base::getDerived());
    }
}

namespace Physica {
    template<Vector V>
    class Traits<ReversalVector<V>> : public Traits<V> {};
}
