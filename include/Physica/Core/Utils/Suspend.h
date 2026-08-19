/*
 * Copyright 2025-2026 Weibo He.
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

    [[nodiscard]] PHYSICA_API std::coroutine_handle<> noop_coroutine() noexcept;
}
