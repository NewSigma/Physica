/*
 * Copyright 2024 WeiBo He.
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

#include <Physica/Utils/Container/Array/Array.h>

namespace Physica::Core {
    template<unsigned int Dim>
    class SiteIndex : public Utils::Array<unsigned char, Dim + 1> {
        using This = SiteIndex<Dim>;
        using Base = Utils::Array<unsigned char, Dim + 1>;
        static_assert(1 <= Dim && Dim <= 3, "[Error]: Invalid Dim");
    public:
        using Base::Base;
        SiteIndex(const This&) = default;
        SiteIndex(This&&) noexcept = default;
        ~SiteIndex() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] unsigned char toIndex1D(This dims) const noexcept;
        using Base::swap;
        /* Getters */
        [[nodiscard]] unsigned char getSite() const noexcept { return (*this)[Dim]; }
        [[nodiscard]] unsigned char getX() const noexcept { return (*this)[0]; }
        [[nodiscard]] unsigned char getY() const noexcept { return (*this)[1]; }
        [[nodiscard]] unsigned char getZ() const noexcept { return (*this)[2]; }
        /* Static members */
        [[nodiscard]] static unsigned char toIndex1D(This dims, This index) noexcept { return index.toIndex1D(dims); }
    };

    template<unsigned int Dim>
    unsigned char SiteIndex<Dim>::toIndex1D(SiteIndex<Dim> dims) const noexcept {
        unsigned char result = (*this)[0];
        for (unsigned int dim = 1; dim <= Dim; ++dim)
            result = result * dims[dim] + (*this)[dim];
        return result;
    }
}
