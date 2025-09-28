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
#include "Physica/Core/Parallel/ThreadPool.h"
#include "TaskBase.h"

namespace Physica {
    template<>
    class Task<Thread> : public TaskBase {
        using This = Task<Thread>;
        using Base = TaskBase;

        struct ThreadAwaiter : public std::suspend_always {
            void await_suspend(std::coroutine_handle<> handle) const noexcept {
                ThreadPool::getInstance().schedule(handle);
            }
        };

        struct Promise {
            std::exception_ptr ex = nullptr;
        public:
            Task get_return_object() noexcept { return std::coroutine_handle<Promise>::from_promise(*this); }
            ThreadAwaiter initial_suspend() noexcept { return {}; }
            std::suspend_always final_suspend() noexcept { return {}; }
            void return_void() noexcept {}
            void unhandled_exception() noexcept { ex = std::current_exception(); }
        };
    public:
        using promise_type = Promise;
        using Range = std::pair<unsigned int, unsigned int>;
    public:
        Task() = default;
        Task(std::coroutine_handle<Promise> handle_) : Base(std::move(handle_)) {}
        Task(const Task&) = delete;
        Task(Task&& obj) noexcept = default;
        ~Task() = default;
        /* Operators */
        Task& operator=(Task obj) noexcept {
            swap(obj);
            return *this;
        }
        /* Operations */
        void wait();
        void swap(Task& __restrict obj) noexcept { Base::swap(obj); }
        /* Static members */
        [[nodiscard]] inline static Range splitJob(size_t num_loop, int part, int i) noexcept;
    };

    inline void Task<Thread>::wait() {
        while (!done()) {
            auto handle = ThreadPool::getInstance().steal();
            if (handle)
                handle.resume();
            else
                std::this_thread::yield();
        }

        std::exception_ptr ex = Base::handle<Promise>().promise().ex;
        if (ex)
            std::rethrow_exception(ex);
    }

    inline auto Task<Thread>::splitJob(size_t num_loop, int part, int i) noexcept -> Range {
        assert(0 <= i && i < part);
        const size_t maxLoopPerCore = (num_loop + part - 1) / part;
        const size_t from = i * maxLoopPerCore;
        const size_t to = std::min(from + maxLoopPerCore, num_loop);
        return std::make_pair(from, to);
    }

    template<ExecutePolicy P>
    [[nodiscard]] Task<Concurrent> parallel_for(std::invocable<size_t> auto fn, size_t num_loop) noexcept requires(P == Thread) {
        assert(num_loop > 0);
        Array<Task<Thread>> tasks(num_loop);
        for (size_t i = 0; i < num_loop; ++i) {
            tasks[i] = [](auto fn, size_t i) noexcept -> Task<Thread> {
                fn(i);
                co_return;
            }(fn, i);
        }

        co_await std::suspend_always{};
        for (auto& task : tasks)
            task.wait();
    }

    template<ExecutePolicy P>
    [[nodiscard]] Task<Concurrent> parallel_for(std::invocable<size_t> auto fn, size_t num_loop, int part) noexcept requires(P == Thread) {
        const bool shouldInferPart = part <= 0 || part > num_loop;
        if (shouldInferPart)
            part = std::min<size_t>(ThreadPool::getInstance().getNumThreads(), num_loop);

        using Range = Task<Thread>::Range;
        Array<Task<Thread>> tasks(part);
        for (int i = 0; i < part; ++i) {
            tasks[i] = [](auto fn, Range range) noexcept -> Task<Thread> {
                for (size_t loop = range.first; loop < range.second; ++loop)
                    fn(loop);
                co_return;
            }(fn, Task<Thread>::splitJob(num_loop, part, i));
        }

        co_await std::suspend_always{};
        for (auto& task : tasks)
            task.wait();
    }
}
