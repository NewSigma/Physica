/*
 * Copyright 2024-2026 Weibo He.
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

#include <format>
#include "Physica/Core/Scalar/Scalar.h"

namespace Physica {
    template<class T>
    class DiffCoro : public T {
        static_assert(ReverseDiff<T>, "[Error]: DiffCoro binds compute graph to reverse diffable objects");
        static_assert(!is_codiff<T>::value, "[Error]: Nested DiffCoro is not allowed");
        static_assert(std::is_object<T>::value, "[Error]: Must save the return by value");
        using This = DiffCoro<T>;
        using Base = T;

        class Promise;
    public:
        using promise_type = Promise;
    private:
        std::coroutine_handle<Promise> handle = nullptr;
    public:
        DiffCoro() = default;
        DiffCoro(ReverseDiff auto&& expr) noexcept requires(!is_codiff<decltype(expr)>::value); // Codiff will delegate to the move constructor
        DiffCoro(const This& other) = delete;
        DiffCoro(This&& other) noexcept;
        [[gnu::nodebug]] ~DiffCoro();
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&& obj) noexcept;
        using Base::operator=;
        /* Operations */
        void reverse_final(auto&&... args) noexcept;
        void swap(This& __restrict obj) noexcept;
    private:
        DiffCoro(Promise& p) noexcept;
        /* Operations */
        [[gnu::nodebug]] void reverse_impl() noexcept;
        [[nodiscard]] static This compute(ReverseDiff auto&& expr) noexcept;
    };

    template<class T>
    class DiffCoro<T>::Promise : public suspend_always {
        using This = Promise;

        T* pObj;
    public:
        Promise() = default;
        Promise(const This&) = delete;
        Promise(This&&) noexcept = delete;
        ~Promise() = default;
        /* Operators */
        This& operator=(const This&) noexcept = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void listen(DiffCoro<T>& node) noexcept;

        auto get_return_object() noexcept;
        static DiffCoro<T> get_return_object_on_allocation_failure() noexcept { unreachable("Expect coro frame is small"); }
        auto initial_suspend() noexcept { return suspend_never{}; }
        void await_transform(auto&&) noexcept = delete("[Error]: Differential coroutines that return non-void values are required to suspend by yielding");
        auto final_suspend() noexcept { return suspend_never{}; }
        auto& yield_value(auto&& arg) noexcept;
        void return_void() noexcept {}
        [[noreturn]] void unhandled_exception() { throw; }

        [[nodiscard]] T& await_resume() const noexcept { return *pObj; }
    };
    //////////////////////////////////////////////////////////////////////////////
    template<>
    class DiffCoro<void> {
        using This = DiffCoro<void>;
        class Promise;
    public:
        using promise_type = Promise;
    private:
        std::coroutine_handle<Promise> handle = nullptr;
    public:
        DiffCoro() = default;
        DiffCoro(const This&) = delete;
        DiffCoro(This&& other) noexcept : handle(other.handle) { other.handle = nullptr; }
        [[gnu::nodebug]] ~DiffCoro() {
            if (handle) {
                assert(!handle.done() && "[Error]: Unexpected resume, this is a bug");
                handle.resume();
                handle = nullptr;
            }
        }
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(This& __restrict obj) noexcept {
            assert(this != &obj && "[Error]: Self swap is likely a bug");
            std::swap(handle, obj.handle);
        }
    private:
        DiffCoro(Promise& p) noexcept : handle(std::coroutine_handle<Promise>::from_promise(p)) {}
    };

    class DiffCoro<void>::Promise {
    public:
        auto get_return_object() noexcept { return DiffCoro<void>(*this); };
        static DiffCoro<void> get_return_object_on_allocation_failure() noexcept { unreachable("Expect coro frame is small"); }
        static auto initial_suspend() noexcept { return suspend_never{}; }
        static auto final_suspend() noexcept { return suspend_never{}; }
        void return_void() noexcept {}
        [[noreturn]] static void unhandled_exception() { throw; }
    };

    template<class Predicate, class Operation, class Functor>
    [[nodiscard("[Warn]: Discarding a coroutine")]] auto co_for(Predicate&& pred, Operation&& op, Functor&& func) noexcept(noexcept(func()));

    template<class Predicate, class FuncTrue, class FuncFalse>
    [[nodiscard("[Warn]: Discarding a coroutine")]] auto co_if(Predicate&& pred, FuncTrue&& t, FuncFalse&& f) noexcept(noexcept(t()) && noexcept(f()));

    template<class Predicate, class Functor>
    [[nodiscard("[Warn]: Discarding a coroutine")]] auto co_if(Predicate&& pred, Functor&& func) noexcept(noexcept(func()));
}

namespace Physica {
    template<class T>
    class Traits<DiffCoro<T>> : public Traits<T> {};
}

namespace std {
    template<class T>
    struct formatter<Physica::DiffCoro<T>, char> {
        constexpr auto parse(std::format_parse_context& ctx) {
            return ctx.begin();
        }

        static auto format(const Physica::DiffCoro<T>& obj, auto& ctx) {
            return formatter<T, char>{}.format(obj, ctx);
        }
    };
}

#include "DiffCoroImpl.h"
