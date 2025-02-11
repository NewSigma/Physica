/*
 * Copyright 2022-2024 Weibo He.
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

#include <functional>
#include "Physica/Core/Parallel/Future/ProcessFuture.h"

namespace Physica {
    class PHYSICA_API SubProcess {
    private:
        std::function<void()> task;
        pid_t pid;
        int nice_incr;
    public:
        SubProcess();
        SubProcess(std::function<void()> task_, int nice_incr_);
        SubProcess(const SubProcess&) = delete;
        SubProcess(SubProcess&& process) noexcept = default;
        ~SubProcess() = default;
        /* Operators */
        SubProcess& operator=(SubProcess process) noexcept;
        /* Operations */
        ProcessFuture execute();
        void swap(SubProcess& __restrict process) noexcept;
        /* Getters */
        [[nodiscard]] pid_t getPid() const noexcept { return pid; }
    };
}
