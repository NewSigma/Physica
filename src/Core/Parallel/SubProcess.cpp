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
#include <iostream>
#include <unistd.h>
#include <signal.h>
#include <sys/wait.h>
#include <sys/prctl.h>
#include "Physica/Core/Parallel/SubProcess.h"

namespace Physica::Core::Parallel {
    SubProcess::SubProcess(SubProcess&& process) noexcept : task(std::move(process.task)), pid(process.pid) {}

    SubProcess& SubProcess::operator=(SubProcess process) noexcept {
        swap(process);
        return *this;
    }

    ProcessFuture SubProcess::execute() {
        pid = fork();
        if (pid == -1) {
            std::cerr << "[Error]: Failed to fork process.\n";
            exit(EXIT_FAILURE);
        }
        else if (pid == 0) {
            prctl(PR_SET_PDEATHSIG, SIGTERM);
            task();
            _exit(EXIT_SUCCESS);
        }
        return ProcessFuture(pid);
    }

    void SubProcess::swap(SubProcess& process) noexcept {
        task.swap(process.task);
        std::swap(pid, process.pid);
    }
}
