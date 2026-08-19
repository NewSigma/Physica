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
#include "Physica/Core/Parallel/Task/Task.h"

using namespace Physica;

Task::Task(Task&& obj) noexcept : h(std::exchange(obj.h, nullptr)) {}

Task::~Task() noexcept {
    if (!empty()) {
        [[maybe_unused]] auto ex = wait(std::nothrow);
        assert(ex == nullptr && "[Error]: Exception escaped");
        h.destroy();
    }
}

auto Task::operator co_await() const noexcept -> suspend_transfer {
    assert(!empty());
    return suspend_transfer{.child = h};
}

void Task::wait() {
    if (auto ex = wait(std::nothrow)) [[unlikely]]
        std::rethrow_exception(ex);
}

auto Task::wait(std::nothrow_t) noexcept -> std::exception_ptr {
    auto& p = h.promise();
    while (!done()) {
        auto waiter = p.waiter();
        assert(waiter != nullptr && "[Error]: Unexpected empty waiter, this is a bug");
        waiter(h);
    }
    return p.exception();
}

void Task::swap(Task& __restrict obj) noexcept {
    assert(this != &obj && "[Error]: Self swap is likely a bug");
    std::swap(h, obj.h);
}

bool Task::done() const noexcept {
    return h.promise().done();
}

bool Task::empty() const noexcept {
    return h == nullptr;
}

Task::promise_type::promise_type() noexcept : continuation(noop_coroutine()) {}

auto Task::promise_type::initial_suspend() noexcept -> suspend_never {
    return {};
}

auto Task::promise_type::final_suspend() noexcept -> suspend_final {
    return {};
}

void Task::promise_type::unhandled_exception() noexcept {
    ex = std::current_exception();
}

auto Task::promise_type::schedule(Handle parent) noexcept -> std::coroutine_handle<> {
    auto expected = noop_coroutine();
    if (continuation.compare_exchange_strong(expected, parent, std::memory_order_acq_rel))
        parent.promise().onWait.store(waiter(), std::memory_order_release);
    else
        assert(expected == nullptr && "[Error]: child has more than one parent, this is a bug");
    return expected;
}

auto Task::promise_type::schedule(std::nullptr_t) noexcept -> std::coroutine_handle<> {
    return continuation.exchange(nullptr, std::memory_order_acq_rel);
}

auto Task::promise_type::waiter() const noexcept -> Callback {
    return onWait.load(std::memory_order_acquire);
}

bool Task::promise_type::done() const noexcept {
    return continuation.load(std::memory_order_acquire) == nullptr;
}

auto Task::promise_type::exception() noexcept -> std::exception_ptr {
    return std::move(ex);
}

auto Task::promise_type::suspend_final::await_suspend(std::coroutine_handle<> self) noexcept -> std::coroutine_handle<> {
    return Handle::from_address(self.address()).promise().schedule(nullptr);
}

bool Task::suspend_transfer::await_ready() const noexcept {
    return child.promise().done();
}

bool Task::suspend_transfer::await_suspend(std::coroutine_handle<> parent) const noexcept {
    assert(parent != nullptr);
    bool suspend = child.promise().schedule(Handle::from_address(parent.address())) != nullptr;
    return suspend;
}

void Task::suspend_transfer::await_resume() const {
    assert(child.promise().done());
    if (auto ex = child.promise().exception()) [[unlikely]]
        std::rethrow_exception(ex);
}
