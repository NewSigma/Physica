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
#include "Physica/Core/Utils/Builtin.h"

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

std::exception_ptr Task::wait(std::nothrow_t) noexcept {
    while (!done())
        __builtin_ia32_pause();
    return h.promise().exception();
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

Task Task::promise_type::get_return_object() noexcept {
    return Task(handle());
}

Task Task::promise_type::get_return_object_on_allocation_failure() noexcept {
    unreachable("Expect coro frame is small");
}

auto Task::promise_type::initial_suspend() noexcept -> suspend_never {
    return {};
}

auto Task::promise_type::final_suspend() noexcept -> suspend_final {
    return {};
}

void Task::promise_type::return_void() noexcept {}

void Task::promise_type::unhandled_exception() noexcept {
    ex = std::current_exception();
}

std::coroutine_handle<> Task::promise_type::schedule(std::coroutine_handle<> todo) noexcept {
    return continuation.exchange(todo, std::memory_order_acq_rel);
}

std::exception_ptr Task::promise_type::exception() noexcept {
    return std::move(ex);
}

bool Task::promise_type::done() const noexcept {
    return continuation.load(std::memory_order_acquire) == nullptr;
}

auto Task::promise_type::handle() noexcept -> Handle {
    return Handle::from_promise(*this);
}

auto Task::promise_type::suspend_final::await_suspend(std::coroutine_handle<> h) noexcept -> std::coroutine_handle<> {
    return Handle::from_address(h.address()).promise().schedule(nullptr);
}

bool Task::suspend_transfer::await_ready() const noexcept {
    return child.promise().done();
}

bool Task::suspend_transfer::await_suspend(std::coroutine_handle<> parent) const noexcept {
    const auto todo = child.promise().schedule(parent);
    bool done = todo == nullptr;
    if (done)
        return false;
    assert(todo == noop_coroutine() && "[Error]: child has more than one parent, this is a bug");
    return true;
}

void Task::suspend_transfer::await_resume() const {
    if (auto ex = child.promise().exception()) [[unlikely]]
        std::rethrow_exception(ex);
}
