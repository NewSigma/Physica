/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Macro.h"
#include "Physica/Core/Utils/Handle.h"

namespace Physica {
    class PHYSICA_API ProcessFuture {
        using This = ProcessFuture;
    public:
        using Handle = Handle<HandleType::PID>;
    private:
        Handle pid;
        int error;
        bool finished;
        bool isValid;
    public:
        ProcessFuture();
        ProcessFuture(Handle pid_);
        ProcessFuture(const This&) = default;
        ProcessFuture(This&&) noexcept = default;
        ~ProcessFuture() = default;
        /* Operators */
        This& operator=(This future) noexcept { swap(future); return *this; }
        /* Operations */
        [[nodiscard]] int wait();
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] Handle getPID() const noexcept { return pid; }
        [[nodiscard]] bool valid() const noexcept { return isValid; }
    };
}
