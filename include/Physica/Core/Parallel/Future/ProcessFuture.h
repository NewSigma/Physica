/*
 * Copyright 2022-2026 Weibo He.
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

#include "Physica/Macro.h"
#include "Physica/Core/Utils/Handle.h"

namespace Physica {
    class PHYSICA_API ProcessFuture {
        using This = ProcessFuture;
    public:
        using Handle = Handle<HandleType::PID>;
    private:
        Handle pid;
        int error;
        bool finished;
        bool isValid;
    public:
        ProcessFuture();
        ProcessFuture(Handle pid_);
        ProcessFuture(const This&) = default;
        ProcessFuture(This&&) noexcept = default;
        ~ProcessFuture() = default;
        /* Operators */
        This& operator=(This future) noexcept { swap(future); return *this; }
        /* Operations */
        [[nodiscard]] int wait();
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] Handle getPID() const noexcept { return pid; }
        [[nodiscard]] bool valid() const noexcept { return isValid; }
    };
}
