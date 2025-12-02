/*
 * Copyright 2024 Weibo He.
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

#include "Physica/Core/Utils/Container/Array.h"

namespace Physica {
    template<int Dim>
    class SiteIndex : public Array<size_t, Dim + 1> {
        using This = SiteIndex<Dim>;
        using Base = Array<size_t, Dim + 1>;
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
        [[nodiscard]] This shift(int shiftDim, size_t delta, size_t period) const noexcept;
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
    auto SiteIndex<Dim>::shift(int shiftDim, size_t delta, size_t period) const noexcept -> This {
        assert(shiftDim < Dim && "[Error]: Invalid dim");
        This result = *this;
        result[shiftDim] = (result[shiftDim] + delta) % period;
        return result;
    }

    template<int Dim>
    size_t SiteIndex<Dim>::toIndex1D(This dims) const noexcept {
        return Base::toIndex1D(dims, *this);
    }
}
