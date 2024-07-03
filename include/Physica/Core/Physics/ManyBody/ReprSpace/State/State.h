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

#include <cassert>
#include <utility>
#include <stdint.h>
#include <Physica/Core/MultiPrecision/BasicImpl/Util/Bitwise.h>
#include <Physica/Core/Physics/ManyBody/SiteIndex.h>

namespace Physica {
    template<class T> class Traits;
}

namespace Physica::Core {
    template<class Derived>
    class State {
    public:
        constexpr static unsigned int Dim = Traits<Derived>::Dim;
        constexpr static unsigned int NumSite = Traits<Derived>::NumSite;
        constexpr static unsigned int SiteDOF = Traits<Derived>::SiteDOF;
        using IndexType = SiteIndex<Dim>;
    public:
        constexpr static size_t calcFullNumState() noexcept;
    };

    template<class Derived>
    constexpr size_t State<Derived>::calcFullNumState() noexcept {
        size_t result = 1;
        for (unsigned int i = 0; i < NumSite; ++i)
            result *= SiteDOF;
        return result;
    }
}
