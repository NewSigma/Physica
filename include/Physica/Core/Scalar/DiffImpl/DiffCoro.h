/*
 * Copyright 2024-2025 Weibo He.
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

#include "Physica/Core/Scalar/Scalar.h"

namespace Physica::Core {
    template<class T> requires(std::is_reference<T>::value)
    using LazyDestroy = std::conditional<std::is_rvalue_reference<T>::value, std::remove_reference_t<T>, T&>::type;

    template<class Base>
    class DiffCoro : public Base {
        static_assert(ReverseDiff<Base>, "[Error]: DiffCoro save compute graph for reverse diffable objects");
        static_assert(!IsCoDiff<Base>::value, "[Error]: Nested DiffCoro is not allowed");
        static_assert(std::is_object<Base>::value, "[Error]: Must save the return by value");
        using This = DiffCoro<Base>;

        class Promise {
            struct suspend_yield : public std::suspend_always {
                std::coroutine_handle<Promise> handle;

                void await_suspend(std::coroutine_handle<Promise> handle_) noexcept {
                    handle = std::move(handle_);
                }

                [[nodiscard]] const Base await_resume() const noexcept {
                    return Base(std::move(handle.promise().obj));
                }
            };
        public:
            Base obj;
        public:
            Promise() = default;
            Promise(const Promise&) = default;
            Promise(Promise&&) noexcept = default;
            ~Promise() = default;
            /* Operators */
            Promise& operator=(Promise obj) noexcept { swap(obj); return *this; }
            /* Operations */
            auto get_return_object() { return std::coroutine_handle<Promise>::from_promise(*this); };
            std::suspend_never initial_suspend() noexcept { return {}; }
            std::suspend_always final_suspend() noexcept { return {}; }
            template<class T>
            auto yield_value(T&& arg) noexcept {
                obj = Base(std::forward<T>(arg));
                return suspend_yield{};
            }
            void return_void() noexcept {}
            void unhandled_exception() { throw std::current_exception(); }

            void swap(This& __restrict obj_) noexcept { obj.swap(obj_); }
        };
    public:
        using promise_type = Promise;
        using typename Base::ScalarType;
    private:
        std::coroutine_handle<Promise> handle = nullptr;
    public:
        DiffCoro() = default;
        DiffCoro(std::coroutine_handle<Promise> handle_) noexcept;
        template<ReverseDiff T>
        DiffCoro(T&& x) noexcept requires(!IsCoDiff<T>::value);
        DiffCoro(const This& other) = delete;
        DiffCoro(This&& other) noexcept;
        ~DiffCoro();
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&& obj) noexcept;
        using Base::operator=;
        /* Operations */
        template<class T>
        void reverse_final(T&& x) noexcept;
        void swap(This& __restrict obj) noexcept;
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

#include "DiffCoroImpl.h"
