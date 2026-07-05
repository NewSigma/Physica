/*
 * Copyright 2021-2026 Weibo He.
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

#include <cassert>
#include <coroutine>
#include <exception>
#include "Physica/Core/Parallel/Parallel.h"

namespace Physica {
    class TaskBase {
        using This = TaskBase;

        std::coroutine_handle<> h = nullptr;
    public:
        TaskBase() = default;
        TaskBase(const This&) = delete;
        ~TaskBase();
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        This& operator=(std::nullptr_t) noexcept;
        [[nodiscard]] inline operator bool() const noexcept;
        /* Operations */
        void resume();
        /* Getters */
        [[nodiscard]] bool empty() const noexcept;
        [[nodiscard]] bool done() const noexcept;
    protected:
        TaskBase(std::coroutine_handle<> h_) : h(h_) {}
        TaskBase(This&& obj) noexcept;
        /* Operations */
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        template<class Promise>
        [[nodiscard]] std::coroutine_handle<Promise> handle() const noexcept;
    };

    inline TaskBase::TaskBase(This&& obj) noexcept : h(std::exchange(obj.h, nullptr)) {}

    inline TaskBase::~TaskBase() {
        if (!empty()) {
            assert(std::uncaught_exceptions() == 0 && "[Error]: Task escapes on unwinding");
            assert(h.done() && "[Error]: Task escapes before finishing");
            h.destroy();
            h = nullptr;
        }
    }

    inline auto TaskBase::operator=(std::nullptr_t) noexcept -> This& {
        h = nullptr;
        return *this;
    }

    inline TaskBase::operator bool() const noexcept {
        return !empty();
    }

    inline void TaskBase::resume() {
        assert(!done());
        h.resume();
    }

    inline bool TaskBase::empty() const noexcept {
        return h == nullptr;
    }

    inline bool TaskBase::done() const noexcept {
        assert(!empty());
        return h.done();
    }

    template<class Promise>
    std::coroutine_handle<Promise> TaskBase::handle() const noexcept {
        return std::coroutine_handle<Promise>::from_address(h.address());
    }

    inline void TaskBase::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(h, obj.h);
    }
}
