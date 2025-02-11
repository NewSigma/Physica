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
#include <utility>
#include "Physica/Core/Utils/SSHWarpper.h"
#include "Physica/Core/Exception/SystemException.h"
#include "Physica/Core/Exception/NoImplException.h"

namespace Physica {
    SSHWarpper::SSHWarpper(std::string hostname_, std::string command_)
            : hostname(std::move(hostname_)), command(std::move(command_)) {
        execute();
    }

    void SSHWarpper::execute() {
    #ifdef __linux__
        int fd[2];
        if (pipe(fd) == -1)
            throw SystemException();

        future = ProcessExecutor::schedule([this, fd]() {
            const int standardErr = dup(STDERR_FILENO);
            close(STDOUT_FILENO);
            close(STDERR_FILENO);
            close(fd[1]);
            if (dup2(fd[0], STDIN_FILENO) != STDIN_FILENO)
                throw SystemException();
            close(fd[0]);
            /* Execute */ {
                execlp("ssh", "ssh", "-tt", hostname.c_str(), nullptr);
            }
            dup2(standardErr, 2);
            perror("[Error]: Failed to execute SSH");
            _exit(EXIT_FAILURE);
        });

        close(fd[0]);
        if (write(fd[1], command.c_str(), command.size()) == -1)
            throw SystemException();
        close(fd[1]);
    #else
        noImpl(__func__);
    #endif
    }

    void SSHWarpper::swap(SSHWarpper& __restrict ssh) noexcept {
        assert(this != &ssh && "[Error]: Self swap is likely a bug");
        hostname.swap(ssh.hostname);
        command.swap(ssh.command);
        future.swap(ssh.future);
    }
}
