/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <coroutine>
#include "MetaProgramming.h"
#include "Physica/Macro.h"

namespace Physica {
    template<class T>
    concept Awaitable = requires(T t, std::coroutine_handle<> handle) {
        { t.await_ready() } -> std::same_as<bool>;
        { t.await_suspend(handle) } -> Either<void, bool, std::coroutine_handle<>>;
        { t.await_resume() };
    };

    template<class T>
    concept Waitable = requires(std::coroutine_handle<> waitee) {
        { T::on_wait(waitee) } -> std::same_as<void>;
    };

    [[nodiscard]] decltype(auto) toAwaiter(auto&& expr) noexcept {
        using Expr = decltype(expr);
        if constexpr (Awaitable<Expr>)
            return std::forward<Expr>(expr);
        else
            return std::forward<Expr>(expr).operator co_await();
    }

    [[nodiscard, gnu::const]] PHYSICA_API std::coroutine_handle<> noop_coroutine() noexcept;
    // Clang implements await_suspend using an intrinsic. We provide the necessary information to help with optimization.
    template<bool Ready>
    struct StaticSuspend {
        // Discard the handle, suspension does not escape coroutine frame
        struct NoEscape {
            NoEscape(auto) {}
        };

        // Static member function, suspension does not escape awaiter
        constexpr static bool await_ready() noexcept { return Ready; }
        constexpr static void await_suspend(NoEscape) noexcept {}
        constexpr static void await_resume() noexcept {}
    };

    using suspend_always = StaticSuspend<false>;
    using suspend_never = StaticSuspend<true>;
}
