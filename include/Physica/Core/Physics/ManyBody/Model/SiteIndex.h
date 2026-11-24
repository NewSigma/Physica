/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/IndexVar/Index.h"

namespace Physica {
    template<int Dim>
    class SiteIndex : public Index<Dim + 1> {
        using This = SiteIndex<Dim>;
        using Base = Index<Dim + 1>;
        static_assert(1 <= Dim && Dim <= 3, "[Error]: Invalid Dim");
    public:
        using Base::Base;
        SiteIndex(const Base& base);
        SiteIndex(const This&) = default;
        SiteIndex(This&&) noexcept = default;
        ~SiteIndex() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] This shift(int shiftDim, ssize_t delta, size_t period) const noexcept;
        [[nodiscard]] size_t toIndex1D(This dims) const noexcept;

        using Base::swap;
        /* Getters */
        [[nodiscard]] size_t getSite() const noexcept { return (*this)[Dim]; }
        [[nodiscard]] size_t getX() const noexcept { return (*this)[0]; }
        [[nodiscard]] size_t getY() const noexcept { return (*this)[1]; }
        [[nodiscard]] size_t getZ() const noexcept { return (*this)[2]; }
        /* Static members */
        using Base::toIndex1D;
    };

    template<int Dim>
    SiteIndex<Dim>::SiteIndex(const Base& base) : Base(base) {}

    template<int Dim>
    auto SiteIndex<Dim>::shift(int shiftDim, ssize_t delta, size_t period) const noexcept -> This {
        assert(shiftDim < Dim && "[Error]: Invalid dim");
        assert(static_cast<size_t>(std::abs(delta)) < period);
        This result = *this;
        result[shiftDim] = (result[shiftDim] + (ssize_t(period) + delta)) % period;
        return result;
    }

    template<int Dim>
    size_t SiteIndex<Dim>::toIndex1D(This dims) const noexcept {
        return Base::toIndex1D(dims, *this);
    }
}
