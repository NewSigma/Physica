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

#include <exception>
#include "Physica/Core/Parallel/ThreadPool.h"
#include "TaskBase.h"

namespace Physica {
    template<>
    class Task<Concurrent> : public TaskBase {
        using This = Task<Concurrent>;
        using Base = TaskBase;

        struct Promise {
            std::exception_ptr ex = nullptr;
        public:
            Task get_return_object() noexcept { return std::coroutine_handle<Promise>::from_promise(*this); }
            static Task get_return_object_on_allocation_failure() noexcept { unreachable("Expect coro frame is small"); }
            auto initial_suspend() noexcept { return suspend_never{}; }
            auto final_suspend() noexcept { return suspend_always{}; }
            void return_void() noexcept {}
            void unhandled_exception() noexcept { ex = std::current_exception(); }
        };
    public:
        using promise_type = Promise;
    public:
        Task() = default;
        Task(std::coroutine_handle<Promise> handle) : Base(handle) {}
        Task(const Task&) = delete;
        Task(Task&& obj) noexcept = default;
        ~Task();
        /* Operators */
        Task& operator=(Task obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] std::exception_ptr wait(std::nothrow_t) noexcept;
        void wait();
        void swap(Task& __restrict obj) noexcept { Base::swap(obj); }
        /* Getters */
        [[nodiscard]] std::exception_ptr getException() const noexcept;
    };

    inline Task<Concurrent>::~Task() {
        if (!empty()) {
            [[maybe_unused]] auto ex = wait(std::nothrow);
            assert(ex == nullptr && "Exception escape");
        }
    }

    inline std::exception_ptr Task<Concurrent>::wait(std::nothrow_t) noexcept {
        while (!done()) {
            resume();
            ThreadPool::spin();
        }
        return getException();
    }

    inline void Task<Concurrent>::wait() {
        if (auto ex = wait(std::nothrow)) [[unlikely]]
            std::rethrow_exception(ex);
    }

    inline std::exception_ptr Task<Concurrent>::getException() const noexcept {
        return handle<Promise>().promise().ex;
    }
}
