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

#include "TaskBase.h"
#include "Physica/Core/Utils/Suspend.h"

namespace Physica {
    template<>
    class Task<Sequential> : public TaskBase {
        using This = Task<Sequential>;
        using Base = TaskBase;

        struct Promise {
            Task get_return_object() noexcept { return {}; }
            auto initial_suspend() noexcept { return suspend_never{}; }
            void await_transform(auto&&) noexcept = delete;
            auto final_suspend() noexcept { return suspend_never{}; }
            void return_void() noexcept {}
            [[noreturn]] void unhandled_exception() { throw; }
        };
    public:
        using promise_type = Promise;
    public:
        Task() = default;
        Task(const Task&) = delete;
        Task(Task&& obj) noexcept = default;
        ~Task() = default;
        /* Operators */
        Task& operator=(Task obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void wait() {}
        void swap(Task& __restrict obj) noexcept { Base::swap(obj); }
        /* Getters */
        [[nodiscard]] constexpr static bool done() noexcept { return true; }
    };
    /**
     * always_inline: Avoid the destructor overhead of \class Task
     * nodebug: This function is less interesting
     */
    template<ExecutePolicy P>
    [[gnu::always_inline, gnu::nodebug]] Task<Sequential> parallel_for(auto fn, size_t num_loop) requires(P == Sequential) {
        assert(num_loop > 0);
        for (size_t i = 0; i < num_loop; ++i)
            fn(i);
        return {};
    }

    template<ExecutePolicy P>
    [[gnu::always_inline, gnu::nodebug]] Task<Sequential> parallel_for(auto fn, size_t num_loop, int) requires(P == Sequential) {
        for (size_t i = 0; i < num_loop; ++i)
            fn(i);
        return {};
    }
}
