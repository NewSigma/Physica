/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
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
        [[nodiscard]] constexpr size_t getLength() const noexcept { return v.getLength(); }
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
