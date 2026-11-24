/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <cstddef>
#include "Physica/Macro.h"

namespace Physica {
    template<class Derived>
    class State {
    public:
        constexpr static int Dim = Traits<Derived>::Dim;
        constexpr static int NumSite = Traits<Derived>::NumSite;
        constexpr static int SiteDOF = Traits<Derived>::SiteDOF;

        static_assert(1 <= Dim && Dim <= 3, "[Error]: Invalid Dim");
        static_assert(NumSite >= 0);
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
