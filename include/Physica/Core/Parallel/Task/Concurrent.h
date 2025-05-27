/*
 * Copyright 2025 Weibo He.
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
            std::suspend_never initial_suspend() noexcept { return {}; }
            std::suspend_always final_suspend() noexcept { return {}; }
            void return_void() noexcept {}
            void unhandled_exception() noexcept { ex = std::current_exception(); }
        };
    public:
        using promise_type = Promise;
    public:
        Task() = default;
        Task(std::coroutine_handle<Promise> handle_) : Base(std::move(handle_)) {}
        Task(const Task&) = delete;
        Task(Task&& obj) noexcept = default;
        ~Task() = default;
        /* Operators */
        Task& operator=(Task obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void wait();
        void swap(Task& __restrict obj) noexcept { Base::swap(obj); }
    };

    inline void Task<Concurrent>::wait() {
        auto handle = Base::handle<Promise>();
        while (!this->done())
            handle.resume();

        std::exception_ptr ex = handle.promise().ex;
        if (ex)
            std::rethrow_exception(ex);
    }
}
