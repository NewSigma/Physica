/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <mutex>
#include "Physica/Macro.h"

namespace Physica::Internal {
    /**
     * Protect thread safe
     */
    class ThreadGuardFFTW final {
        using This = ThreadGuardFFTW;
    public:
        std::mutex globalMutex;
    public:
        ThreadGuardFFTW(const This&) = delete;
        ThreadGuardFFTW(This&&) noexcept = delete;
        ~ThreadGuardFFTW() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Static members */
        [[nodiscard]] PHYSICA_API static This& getInstance() noexcept;
    private:
        constexpr ThreadGuardFFTW() = default;
    };
}
