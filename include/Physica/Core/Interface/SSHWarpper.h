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

#include "Physica/Core/Parallel/Executor/ProcessExecutor.h"

namespace Physica::Core {
    class SSHWarpper {
        std::string hostname;
        std::string command;
        mutable Parallel::ProcessFuture future;
    public:
        SSHWarpper() = default;
        SSHWarpper(std::string hostname_, std::string command_);
        SSHWarpper(const SSHWarpper&) = default;
        SSHWarpper(SSHWarpper&&) noexcept = default;
        ~SSHWarpper() = default;
        /* Operators */
        SSHWarpper& operator=(SSHWarpper ssh) noexcept;
        /* Operations */
        void execute();
        /* Getters */
        Parallel::ProcessFuture& getFuture() noexcept { return future; }
        /* Helpers */
        void swap(SSHWarpper& ssh) noexcept;
    };
}
