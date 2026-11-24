/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Parallel/Parallel.h"
#include <cassert>
#include <cstdint>
#include <exception>
#include <concepts>
#include <new>

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
    /**
     * FIXME: GPU fallback is deprecated and would be removed in the future
     */
    template<ExecutePolicy P>
    [[gnu::always_inline, gnu::nodebug]] EmptyTask schedule(std::invocable<> auto fn) requires(P == Sequential || P == GPU) {
        fn();
        return {};
    }

    template<ExecutePolicy P>
    [[gnu::always_inline, gnu::nodebug]] EmptyTask parallel_for(auto fn, size_t num_loop) requires(P == Sequential || P == GPU) {
        for (size_t i = 0; i < num_loop; ++i)
            fn(i);
        return {};
    }

    template<ExecutePolicy P>
    [[gnu::always_inline, gnu::nodebug]] EmptyTask parallel_for(auto fn, size_t num_loop, size_t) requires(P == Sequential || P == GPU) {
        for (size_t i = 0; i < num_loop; ++i)
            fn(i);
        return {};
    }
}
