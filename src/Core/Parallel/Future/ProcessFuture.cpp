/*
 * Copyright 2022-2025 Weibo He.
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
#include <cassert>
#include <signal.h>
#ifdef __linux__
    #include <sys/wait.h>
 #endif
#include "Physica/Core/Parallel/Future/ProcessFuture.h"
#include "Physica/Core/Exception/SystemException.h"
#include "Physica/Core/Exception/NoImplException.h"

using namespace Physica;

ProcessFuture::ProcessFuture()
        : error(-1), isValid(false) {}

ProcessFuture::ProcessFuture(pid_t pid_)
        : pid(pid_), error(-1), finished(false), isValid(true) {}

int ProcessFuture::wait() {
    if (!isValid)
        throw std::future_error(std::future_errc::no_state);

    if (finished)
        return error;
#ifdef __linux__
    int status;
    pid_t endPid = waitpid(pid, &status, 0);
    if (endPid <= 0) {
        fprintf(stderr, "[Error]: Failed to wait for chile processes.\n");
        throw SystemException();
    }
    finished = true;

    if (WIFEXITED(status))
        error = WEXITSTATUS(status);
#else
    noImpl(__func__);
#endif
    return error;
}

void ProcessFuture::swap(ProcessFuture& __restrict obj) noexcept {
    assert(this != &obj && "[Error]: Self swap is likely a bug");
    std::swap(pid, obj.pid);
    std::swap(error, obj.error);
    std::swap(finished, obj.finished);
    std::swap(isValid, obj.isValid);
}
