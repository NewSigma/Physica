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
#include "Physica/Core/Exception/SyscallException.h"

namespace Physica::Core {
    SubProcess::SubProcess() : task(), pid(-1), nice_incr(0) {}

    SubProcess::SubProcess(std::function<void()> task_, int nice_incr_) : task(std::move(task_)), pid(-1), nice_incr(nice_incr_) {}

    SubProcess& SubProcess::operator=(SubProcess process) noexcept {
        swap(process);
        return *this;
    }

    ProcessFuture SubProcess::execute() {
        pid = fork();
        if (pid == -1)
            throw SyscallException();
        else if (pid == 0) {
            prctl(PR_SET_PDEATHSIG, SIGTERM);
            /* Set nice */ {
                errno = 0;
                const int code = nice(nice_incr);
                if (code == -1 && errno != 0)
                    perror("[Warning]: Failed to process priority");
            }
            task();
            _exit(EXIT_SUCCESS);
        }
        return ProcessFuture(pid);
    }

    void SubProcess::swap(SubProcess& process) noexcept {
        task.swap(process.task);
        std::swap(pid, process.pid);
        std::swap(nice_incr, process.nice_incr);
    }
}
