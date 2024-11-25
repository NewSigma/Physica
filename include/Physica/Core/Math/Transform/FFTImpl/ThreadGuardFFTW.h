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

#include <mutex>
#include "Physica/Macro.h"

namespace Physica::Core::Internal {
    /**
     * Protect thread safe
     */
    class ThreadGuardFFTW final {
        using This = ThreadGuardFFTW;
    public:
        std::mutex globalMutex;
    public:
        ThreadGuardFFTW(const This&) = delete;
        ThreadGuardFFTW(This&&) noexcept = delete;
        ~ThreadGuardFFTW() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Static members */
        [[nodiscard]] PHYSICA_API static This& getInstance() noexcept;
    private:
        constexpr ThreadGuardFFTW() = default;
    };
}
