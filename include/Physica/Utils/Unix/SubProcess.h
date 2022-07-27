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

#include "UnixHelper.h"
#include <functional>

namespace Physica::Utils {
    class SubProcess {
    private:
        std::function<void()> task;
        pid_t pid;
    public:
        SubProcess() : task(), pid(-1) {}
        SubProcess(std::function<void()> _task) : task(std::move(_task)), pid(-1) {}
        SubProcess(const SubProcess&) = delete;
        SubProcess(SubProcess&& process) noexcept;
        ~SubProcess() = default;
        /* Operators */
        SubProcess& operator=(SubProcess process) noexcept;
        /* Operations */
        void execute();
        pid_t wait(const char* errorMsg);
        /* Getters */
        [[nodiscard]] pid_t getPid() const noexcept { return pid; }
        /* Helpers */
        void swap(SubProcess& process) noexcept;
    };
}
