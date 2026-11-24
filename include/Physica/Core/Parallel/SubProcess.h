/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <functional>
#include "Physica/Core/Parallel/Future/ProcessFuture.h"

namespace Physica {
    class PHYSICA_API SubProcess {
    public:
        using Handle = Handle<HandleType::PID>;
    private:
        std::function<void()> task;
        Handle pid;
        int nice_incr;
    public:
        SubProcess();
        SubProcess(std::function<void()> task_, int nice_incr_);
        SubProcess(const SubProcess&) = delete;
        SubProcess(SubProcess&& process) noexcept = default;
        ~SubProcess() = default;
        /* Operators */
        SubProcess& operator=(SubProcess process) noexcept;
        /* Operations */
        ProcessFuture execute();
        void swap(SubProcess& __restrict process) noexcept;
        /* Getters */
        [[nodiscard]] Handle getPid() const noexcept { return pid; }
    };
}
