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

#include "Physica/Utils/Template/CRTPBase.h"

namespace Physica {
    template<class T> class Traits;
}

namespace Physica::Core {
    template<class Derived>
    class ReprBasis : public Utils::CRTPBase<Derived> {
        using Base = Utils::CRTPBase<Derived>;
    public:
        using StateType = typename Traits<Derived>::StateType;
        constexpr static unsigned int Dim = Traits<Derived>::Dim;
    public:
        ~ReprBasis() = default;
        /* Operators */
        [[nodiscard]] StateType operator[](size_t index) const noexcept { return Base::getDerived()[index]; }
        [[nodiscard]] size_t operator[](StateType state) const noexcept { return Base::getDerived()[state]; }
        /* Getters */
        [[nodiscard]] size_t getNumState() const noexcept { return Base::getDerived().getSize(); }
    protected:
        ReprBasis() = default;
        ReprBasis(const ReprBasis&) = default;
        ReprBasis(ReprBasis&&) noexcept = default;
        /* Operators */
        ReprBasis& operator=(const ReprBasis&) = default;
        ReprBasis& operator=(ReprBasis&&) noexcept = default;
    };
}
