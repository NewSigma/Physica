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

#include "Physica/Core/Parallel/Parallel.h"
#include <cassert>
#include <cstdint>
#include <exception>
#include <concepts>

namespace Physica {
    /**
     * \class EmptyTask mimics the interface of \class Task with no state; all operations regress to no-op.
     */
    class EmptyTask final {
        using This = EmptyTask;
    public:
        constexpr EmptyTask() = default;
        constexpr EmptyTask(const This&) = delete;
        constexpr EmptyTask(This&&) noexcept = default;
        constexpr ~EmptyTask() = default;
        /* Operators */
        constexpr This& operator=(const This&) = delete;
        constexpr This& operator=(This&&) noexcept = default;
        /* Operations */
        constexpr void wait() const noexcept {}
        [[nodiscard]] constexpr static std::exception_ptr wait(std::nothrow_t) noexcept { return nullptr; }

        constexpr void swap(This&) noexcept {}
        /* Getters */
        [[nodiscard]] constexpr static bool empty() noexcept { return true; }
        [[nodiscard]] constexpr static bool done() noexcept { return true; }
    };
    
    template<ExecutePolicy P>
    [[gnu::always_inline, gnu::nodebug]] EmptyTask schedule(std::invocable<> auto fn) requires(P == Sequential) {
        fn();
        return {};
    }

    template<ExecutePolicy P>
    [[gnu::always_inline, gnu::nodebug]] EmptyTask parallel_for(auto fn, size_t num_loop) requires(P == Sequential) {
        for (size_t i = 0; i < num_loop; ++i)
            fn(i);
        return {};
    }

    template<ExecutePolicy P>
    [[gnu::always_inline, gnu::nodebug]] EmptyTask parallel_for(auto fn, size_t num_loop, size_t) requires(P == Sequential) {
        for (size_t i = 0; i < num_loop; ++i)
            fn(i);
        return {};
    }
}
