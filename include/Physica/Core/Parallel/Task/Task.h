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
#include "Physica/Core/Utils/Builtin.h"
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
        using Callback = void (*)(std::coroutine_handle<>);

        struct suspend_final;

        std::atomic<std::coroutine_handle<>> continuation;
        std::atomic<Callback> onWait = nullptr;
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
        [[nodiscard]] Task get_return_object() noexcept { return Task(Handle::from_promise(*this)); }
        [[nodiscard]] static Task get_return_object_on_allocation_failure() noexcept { unreachable("Expect coro frame is small"); }
        [[nodiscard]] static suspend_never initial_suspend() noexcept;
        [[nodiscard]] static suspend_final final_suspend() noexcept;
        [[nodiscard]] auto await_transform(auto&& expr) noexcept;
        static void return_void() noexcept {}
        void unhandled_exception() noexcept;

        [[nodiscard]] auto schedule(Handle parent) noexcept -> std::coroutine_handle<>;
        [[nodiscard]] auto schedule(std::nullptr_t) noexcept -> std::coroutine_handle<>;
        /* Getters */
        [[nodiscard]] auto waiter() const noexcept -> Callback;
        [[nodiscard]] bool done() const noexcept;
        [[nodiscard]] auto exception() noexcept -> std::exception_ptr;
    };

    struct Task::promise_type::suspend_final : public suspend_always {
        [[nodiscard]] static auto await_suspend(std::coroutine_handle<> self) noexcept -> std::coroutine_handle<>;
    };

    auto Task::promise_type::await_transform(auto&& expr) noexcept {
        auto awaiter = toAwaiter(std::forward<decltype(expr)>(expr));
        struct WaitWrapper {
            using Awaiter = decltype(awaiter);
            Awaiter awaiter;
            This* promise;

            bool await_ready() {
                return awaiter.await_ready();
            }

            auto await_suspend(std::coroutine_handle<> h) {
                if constexpr (Waitable<Awaiter>)
                    promise->onWait.store(&Awaiter::on_wait, std::memory_order_release);
                else
                    static_assert(std::same_as<suspend_transfer, Awaiter>, "[Error]: Must co_await a Task or Waitable, otherwise coroutine leaks");
                return awaiter.await_suspend(h);
            }

            decltype(auto) await_resume() {
                return awaiter.await_resume();
            }
        };
        return WaitWrapper{.awaiter = std::move(awaiter), .promise = this};
    }
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
