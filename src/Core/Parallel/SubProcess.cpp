/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#ifdef __linux__

#include <cassert>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <csignal>
#include <sys/wait.h>
#include <sys/prctl.h>
#include <sys/types.h>
#include "Physica/Core/Parallel/SubProcess.h"
#include "Physica/Core/Exception/Syscall.h"

using namespace Physica;

SubProcess::SubProcess()
        : pid(Handle(-1)), nice_incr(0) {}

SubProcess::SubProcess(std::function<void()> task_, int nice_incr_)
        : task(std::move(task_)), pid(Handle(-1)), nice_incr(nice_incr_) {}

SubProcess& SubProcess::operator=(SubProcess process) noexcept {
    swap(process);
    return *this;
}

ProcessFuture SubProcess::execute() {
    pid = Handle(fork());
    check_syscall(pid_t(pid) == -1);
    if (pid_t(pid) == 0) {
        prctl(PR_SET_PDEATHSIG, SIGTERM);
        /* Set nice */ {
            errno = 0;
            const int code = nice(nice_incr);
            if (code == -1 && errno != 0)
                perror("[Warn]: Failed to process priority");
        }
        task();
        _exit(EXIT_SUCCESS);
    }
    return ProcessFuture(pid);
}

void SubProcess::swap(SubProcess& __restrict process) noexcept {
    assert(this != &process && "[Error]: Self swap is likely a bug");
    task.swap(process.task);
    pid.swap(process.pid);
    std::swap(nice_incr, process.nice_incr);
}

#endif
