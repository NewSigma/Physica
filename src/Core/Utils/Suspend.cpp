/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Utils/Suspend.h"

namespace {
    const std::coroutine_handle<> NoopCoro = []() noexcept {
    #ifdef __clang__
        return std::coroutine_handle<>::from_address(__builtin_coro_noop());
    #else
        return std::noop_coroutine();
    #endif
    }();
}

namespace Physica {
    std::coroutine_handle<> noop_coroutine() noexcept {
        return NoopCoro;
    }
}
