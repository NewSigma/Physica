/*
 * Copyright 2022 WeiBo He.
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
#include "FutureGroup.h"

namespace Physica::Core {
    class ProcessFuture {
        pid_t pid;
        bool finished;
        bool isValid;
    public:
        ProcessFuture();
        ProcessFuture(pid_t pid_);
        ProcessFuture(const ProcessFuture&) = default;
        ProcessFuture(ProcessFuture&&) noexcept = default;
        ~ProcessFuture() = default;
        /* Operators */
        ProcessFuture& operator=(ProcessFuture future) noexcept;
        /* Operations */
        void wait(const char* errorMsg);
        /* Getters */
        [[nodiscard]] pid_t getPID() const noexcept { return pid; }
        [[nodiscard]] bool valid() const noexcept { return isValid; }
        /* Helpers */
        void swap(ProcessFuture& __restrict future) noexcept;

        friend class Test;
    };
}
