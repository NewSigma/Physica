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
#include <stdexcept>
#include <future>
#include <unistd.h>
#include <signal.h>
#include <sys/wait.h>
#include "Physica/Core/Parallel/Future/ProcessFuture.h"
#include "Physica/Core/Exception/SyscallException.h"

namespace Physica::Core {
    ProcessFuture::ProcessFuture() : isValid(false) {}

    ProcessFuture::ProcessFuture(pid_t pid_) : pid(pid_), finished(false), isValid(true) {}

    ProcessFuture& ProcessFuture::operator=(ProcessFuture future) noexcept {
        swap(future);
        return *this;
    }

    void ProcessFuture::wait(const char* errorMsg) {
        if (!isValid)
            throw std::future_error(std::future_errc::no_state);

        if (finished)
            return;

        int status;
        pid_t endPid = waitpid(pid, &status, 0);
        if (endPid <= 0) {
            fprintf(stderr, "[Error]: Failed to wait for chile processes.\n");
            throw SyscallException();
        }
        finished = true;

        int error = -1;
        if (WIFEXITED(status))
            error = WEXITSTATUS(status);
        if (error != 0)
            throw std::runtime_error(errorMsg);
    }

    void ProcessFuture::swap(ProcessFuture& future) noexcept {
        std::swap(pid, future.pid);
        std::swap(finished, future.finished);
        std::swap(isValid, future.isValid);
    }
}
