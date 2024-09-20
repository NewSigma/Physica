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

#include <cstddef>
#include <Physica/Macro.h>

namespace Physica::Core {
    template<class Derived>
    class State {
    public:
        constexpr static int Dim = Traits<Derived>::Dim;
        constexpr static int NumSite = Traits<Derived>::NumSite;
        constexpr static int SiteDOF = Traits<Derived>::SiteDOF;

        static_assert(1 <= Dim && Dim <= 3, "[Error]: Invalid Dim");
        static_assert(NumSite > 0);
        static_assert(SiteDOF > 0);
    public:
        constexpr static size_t calcFullNumState() noexcept;
    };

    template<class Derived>
    constexpr size_t State<Derived>::calcFullNumState() noexcept {
        size_t result = 1;
        for (int i = 0; i < NumSite; ++i)
            result *= SiteDOF;
        return result;
    }
}
