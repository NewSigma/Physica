/*
 * Copyright 2022 Weibo He.
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

#include <sys/types.h>
#ifdef _MSC_VER
    //#include <io.h>
    #include <process.h>
    using pid_t = int;
#else
    #include <unistd.h>
#endif
#include "Physica/Macro.h"

namespace Physica {
    class PHYSICA_API ProcessFuture {
        pid_t pid;
        int error;
        bool finished;
        bool isValid;
    public:
        ProcessFuture();
        ProcessFuture(pid_t pid_);
        ProcessFuture(const ProcessFuture&) = default;
        ProcessFuture(ProcessFuture&&) noexcept = default;
        ~ProcessFuture() = default;
        /* Operators */
        ProcessFuture& operator=(ProcessFuture future) noexcept { swap(future); return *this; }
        /* Operations */
        [[nodiscard]] int wait();
        void swap(ProcessFuture& __restrict future) noexcept;
        /* Getters */
        [[nodiscard]] pid_t getPID() const noexcept { return pid; }
        [[nodiscard]] bool valid() const noexcept { return isValid; }
    };
}
