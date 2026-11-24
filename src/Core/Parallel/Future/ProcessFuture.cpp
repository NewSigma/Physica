/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <future>
#include <cassert>
#include <csignal>
#ifdef __linux__
    #include <sys/types.h>
    #include <sys/wait.h>
 #endif
#include "Physica/Core/Parallel/Future/ProcessFuture.h"
#include "Physica/Core/Exception/Syscall.h"

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
    check_syscall(endPid <= 0); // Failed to wait for chile processes
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
