/*
 * Copyright 2021-2025 Weibo He.
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

#include "ThreadPool.h"

namespace Physica {
    enum ParallelPolicy {
        Concurrent,
        Thread
    };

    template<ParallelPolicy P = Concurrent>
    class [[nodiscard]] Task;
    /**
     * \class Task maintains lifetime of async tasks
     */
    template<ParallelPolicy P>
    class Task {
        using This = Task<P>;

        struct ThreadAwaiter : public std::suspend_always {
            void await_suspend(std::coroutine_handle<> handle) const noexcept {
                ThreadPool::getInstance().schedule(handle);
            }
        };

        using InitialAwaiter = std::conditional<P == Concurrent, std::suspend_never, ThreadAwaiter>::type;

        struct Promise {
            std::exception_ptr ex = nullptr;
        public:
            Task get_return_object() noexcept { return std::coroutine_handle<Promise>::from_promise(*this); }
            InitialAwaiter initial_suspend() noexcept { return {}; }
            std::suspend_always final_suspend() noexcept { return {}; }
            void return_void() noexcept {}
            void unhandled_exception() noexcept { ex = std::current_exception(); }
        };
    public:
        using promise_type = Promise;
    private:
        std::coroutine_handle<Promise> handle = nullptr;
    public:
        Task() = default;
        Task(std::coroutine_handle<Promise> handle_) : handle(std::move(handle_)) {}
        Task(const Task&) = delete;
        Task(Task&& obj) noexcept;
        ~Task();
        /* Operators */
        Task& operator=(Task obj) noexcept { swap(obj); return *this; }
        Task& operator=(std::nullptr_t) noexcept;
        [[nodiscard]] inline operator bool() const noexcept;
        /* Operations */
        void wait();
        void swap(Task& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] void* address() const noexcept { return handle.address(); }
        [[nodiscard]] bool done() const noexcept;
    };

    template<ParallelPolicy P>
    inline Task<P>::Task(Task&& obj) noexcept : handle(std::exchange(obj.handle, nullptr)) {}

    template<ParallelPolicy P>
    inline Task<P>::~Task() {
        if (handle) {
            handle.destroy();
            handle = nullptr;
        }
    }

    template<ParallelPolicy P>
    inline auto Task<P>::operator=(std::nullptr_t) noexcept -> This& {
        handle = nullptr;
        return *this;
    }

    template<ParallelPolicy P>
    inline Task<P>::operator bool() const noexcept {
        return bool(handle);
    }

    template<ParallelPolicy P>
    inline void Task<P>::wait() {
        while (!done()) {
            if constexpr (P == Concurrent)
                handle.resume();
            else {
                auto handle = ThreadPool::getInstance().steal();
                if (handle) {
                    handle.resume();
                    if (!handle.done())
                        ThreadPool::getInstance().schedule(handle);
                }
                else
                    std::this_thread::yield();
            }
        }

        std::exception_ptr ex = handle.promise().ex;
        handle.destroy();
        handle = nullptr;
        if (ex)
            std::rethrow_exception(ex);
    }

    template<ParallelPolicy P>
    inline bool Task<P>::done() const noexcept {
        assert(handle != nullptr);
        return handle.done();
    }

    template<ParallelPolicy P>
    inline void Task<P>::swap(Task& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(handle, obj.handle);
    }
}
