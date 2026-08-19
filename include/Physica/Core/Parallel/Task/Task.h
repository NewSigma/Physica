/*
 * Copyright 2026 Weibo He.
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

#include <atomic>
#include <exception>
#include "Physica/Core/Utils/Suspend.h"
#include "Physica/Macro.h"

namespace Physica {
    /**
     * \class Task maintains lifetime of async tasks
     */
    class PHYSICA_API Task final {
    public:
        class promise_type;
    private:
        struct suspend_transfer;
        using Handle = std::coroutine_handle<promise_type>;

        Handle h;
    public:
        Task(const Task&) = delete;
        Task(Task&& obj) noexcept;
        ~Task() noexcept;
        /* Operators */
        Task& operator=(const Task&) = delete;
        Task& operator=(Task&&) noexcept = delete;
        [[nodiscard]] suspend_transfer operator co_await() const noexcept;
        /* Operations */
        void wait();
        [[nodiscard]] std::exception_ptr wait(std::nothrow_t) noexcept;

        void swap(Task& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] bool done() const noexcept;
        [[nodiscard]] bool empty() const noexcept;
    private:
        explicit Task(std::coroutine_handle<promise_type> h_) : h(h_) {}
    };

    class Task::promise_type {
        using This = promise_type;

        struct suspend_final;

        std::atomic<std::coroutine_handle<>> continuation;
        std::exception_ptr ex = nullptr;
    public:
        promise_type() noexcept;
        promise_type(const This&) = delete;
        promise_type(Task&& obj) noexcept = delete;
        ~promise_type() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] Task get_return_object() noexcept;
        [[nodiscard]] static Task get_return_object_on_allocation_failure() noexcept;
        [[nodiscard]] static suspend_never initial_suspend() noexcept;
        [[nodiscard]] static suspend_final final_suspend() noexcept;
        static void return_void() noexcept;
        void unhandled_exception() noexcept;

        [[nodiscard]] std::coroutine_handle<> schedule(std::coroutine_handle<> todo) noexcept;
        /* Getters */
        [[nodiscard]] std::exception_ptr exception() noexcept;
        [[nodiscard]] bool done() const noexcept;
        [[nodiscard]] Handle handle() noexcept;
    };

    struct Task::promise_type::suspend_final : public suspend_always {
        [[nodiscard]] static auto await_suspend(std::coroutine_handle<> h) noexcept -> std::coroutine_handle<>;
    };
    /**
     * Parent task co_await child task, parent is suspended and will be resumed by child
     */
    struct Task::suspend_transfer : public suspend_always {
        Handle child;

        [[nodiscard]] bool await_ready() const noexcept;
        [[nodiscard]] bool await_suspend(std::coroutine_handle<> parent) const noexcept;
        void await_resume() const;
    };
}
