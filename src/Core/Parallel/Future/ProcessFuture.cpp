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
#include <print>
#include <future>
#include <cassert>
#include <csignal>
#ifdef __linux__
    #include <sys/types.h>
    #include <sys/wait.h>
 #endif
#include "Physica/Core/Parallel/Future/ProcessFuture.h"
#include "Physica/Core/Exception/SystemException.h"
#include "Physica/Core/Utils/NoImpl.h"

using namespace Physica;

ProcessFuture::ProcessFuture()
        : pid(-1), error(-1), finished(false), isValid(false) {}

ProcessFuture::ProcessFuture(Handle pid_)
        : pid(pid_), error(-1), finished(false), isValid(true) {}

int ProcessFuture::wait() {
    if (!isValid)
        throw std::future_error(std::future_errc::no_state);

    if (finished)
        return error;
#ifdef __linux__
    int status{};
    const auto endPid = waitpid(pid_t(pid), &status, 0);
    if (endPid <= 0) {
        std::println(stderr, "[Error]: Failed to wait for chile processes.");
        throw SystemException();
    }
    finished = true;

    if (WIFEXITED(status))
        error = WEXITSTATUS(status);
#else
    noImpl();
#endif
    return error;
}

void ProcessFuture::swap(This& __restrict obj) noexcept {
    assert(this != &obj && "[Error]: Self swap is likely a bug");
    pid.swap(obj.pid);
    std::swap(error, obj.error);
    std::swap(finished, obj.finished);
    std::swap(isValid, obj.isValid);
}
