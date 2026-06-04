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

#include <string>
#include "Physica/Core/Parallel/Executor/ProcessExecutor.h"

namespace Physica {
    class PHYSICA_API SSHWrapper {
        using This = SSHWrapper;

        std::string hostname;
        std::string command;
        mutable ProcessFuture future;
    public:
        SSHWrapper() = default;
        SSHWrapper(std::string hostname_, std::string command_);
        SSHWrapper(const This&) = default;
        SSHWrapper(This&&) noexcept = default;
        ~SSHWrapper() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void execute();
        void swap(This& __restrict ssh) noexcept;
        /* Getters */
        ProcessFuture& getFuture() noexcept { return future; }
    };
}
