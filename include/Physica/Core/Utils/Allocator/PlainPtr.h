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

#include <cstdlib>

namespace Physica::Core {
    class PlainPtr {
        using This = PlainPtr;

        void* p = nullptr;
    public:
        PlainPtr() = default;
        explicit PlainPtr(void* p_) : p(p_) {}
        PlainPtr(const This&) = delete;
        PlainPtr(This&&) noexcept = default;
        ~PlainPtr() { std::free(p); }
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = default;
        This& operator=(std::nullptr_t) noexcept { p = nullptr; return *this; }
        /* Getters */
        [[nodiscard]] void* get() const noexcept { return p; }
        [[nodiscard]] inline void* release() noexcept;
    };

    inline void* PlainPtr::release() noexcept {
        void* result = p;
        p = nullptr;
        return result;
    }
}
