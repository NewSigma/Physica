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
#include "Physica/Utils/Unix/SubProcess.h"

namespace Physica::Utils {
    SubProcess::SubProcess(SubProcess&& process) noexcept : task(std::move(process.task)), pid(process.pid) {}

    SubProcess& SubProcess::operator=(SubProcess process) noexcept {
        swap(process);
        return *this;
    }

    void SubProcess::execute() {
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
    }

    pid_t SubProcess::wait(const char* errorMsg) {
        int status;
        pid_t endPid = waitpid(pid, &status, 0);
        pid = -1;
        if (endPid <= 0) {
            fprintf(stderr, "[Error]: Failed to wait for chile processes.\n");
            exit(EXIT_FAILURE);
        }

        int error = -1;
        if (WIFEXITED(status))
            error = WEXITSTATUS(status);
        if (error != 0) {
            fprintf(stderr, "%s\n", errorMsg);
            exit(EXIT_FAILURE);
        }
        return endPid;
    }

    void SubProcess::swap(SubProcess& process) noexcept {
        task.swap(process.task);
        std::swap(pid, process.pid);
    }
}
