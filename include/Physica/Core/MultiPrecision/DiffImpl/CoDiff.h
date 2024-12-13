/*
 * Copyright 2024 Weibo He.
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

#include <coroutine>
#include <exception>
#include <forward_list>
#include "Physica/Core/MultiPrecision/Scalar.h"

namespace Physica::Core {
    template<class T>
    struct IsDiffNode {
        constexpr static bool value = false;
    };

    template<class T>
    struct IsDiffNode<DiffNode<T>> {
        constexpr static bool value = false;
    };

    template<class T> requires(std::is_reference<T>::value)
    using LazyReverse = std::conditional<std::is_rvalue_reference<T>::value, std::remove_reference_t<T>, T&>::type;

    template<class Base>
    class DiffNode : public Base {
        static_assert(ReverseDiff<Base>, "[Error]: DiffNode save compute graph for reverse diffable objects");
        static_assert(!IsDiffNode<Base>::value, "[Error]: Nested DiffNode is not allowed");
        static_assert(std::is_object<Base>::value, "[Error]: Must save the return by value");
        using This = DiffNode<Base>;

        class Impl {
            struct suspend_yield : public std::suspend_always {
                std::coroutine_handle<Impl> handle;

                void await_suspend(std::coroutine_handle<Impl> handle_) noexcept {
                    handle = std::move(handle_);
                }

                [[nodiscard]] Base await_resume() const noexcept {
                    return Base(std::move(handle.promise().obj));
                }
            };
        public:
            Base obj;
        public:
            Impl() = default;
            Impl(const Impl&) = default;
            Impl(Impl&&) noexcept = default;
            ~Impl() = default;
            /* Operators */
            Impl& operator=(Impl obj) noexcept { swap(obj); return *this; }
            /* Operations */
            auto get_return_object() { return std::coroutine_handle<Impl>::from_promise(*this); };
            std::suspend_never initial_suspend() noexcept { return {}; }
            std::suspend_always final_suspend() noexcept { return {}; }
            auto yield_value(Base obj_) noexcept {
                obj = std::move(obj_);
                return suspend_yield{};
            }
            void return_void() noexcept {}
            void unhandled_exception() { throw std::current_exception(); }

            void swap(This& __restrict obj_) noexcept { obj.swap(obj_); }
        };
    public:
        using promise_type = Impl;
        using typename Base::ScalarType;
    private:
        std::coroutine_handle<Impl> handle = nullptr;
    public:
        DiffNode() = default;
        DiffNode(std::coroutine_handle<Impl> handle_);
        DiffNode(const This&) = delete;
        DiffNode(This&& other) noexcept;
        ~DiffNode();
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(This& __restrict obj) noexcept;
    };

    template<class Base>
    DiffNode<Base>::DiffNode(std::coroutine_handle<Impl> handle_) : Base(std::move(handle_.promise().obj)), handle(handle_) {}

    template<class Base>
    DiffNode<Base>::DiffNode(This&& other) noexcept : Base(std::move(other)), handle(other.handle) {
        other.handle = nullptr;
    }

    template<class Base>
    DiffNode<Base>::~DiffNode() {
        if (handle) {
            assert(!handle.done() && "[Error]: Unexpected resume, this is a bug");
            handle.promise().obj = std::move(*this);
            handle.resume();
            assert(handle.done() && "[Error]: Invalid reverse diff");
            handle.destroy();
            handle = nullptr;
        }
    }

    template<class Base>
    void DiffNode<Base>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        std::swap(handle, obj.handle);
    }

    template<class Predicate, class Operation, class Functor>
    [[nodiscard("[Warn]: Discarding a coroutine")]] auto co_for(Predicate&& pred, Operation&& op, Functor&& func) noexcept(func()) {
        static_assert(std::predicate<Predicate>, "[Error]: Invalid predicate");
        using Handle = std::invoke_result<Functor>::type;
        std::forward_list<Handle> result{};
        for (; pred(); op())
            result.push_front(func());
        return result;
    }

    template<class Predicate, class FuncTrue, class FuncFalse>
    [[nodiscard("[Warn]: Discarding a coroutine")]] auto co_if(Predicate&& pred, FuncTrue&& t, FuncFalse&& f) noexcept(noexcept(t()) && noexcept(f())) {
        using Handle1 = std::invoke_result<FuncTrue>::type;
        using Handle2 = std::invoke_result<FuncFalse>::type;
        static_assert(std::is_same<Handle1, Handle2>::value, "[Error]: Inconsistend return value");
        static_assert(std::predicate<Predicate>, "[Error]: Invalid predicate");

        if (pred())
            return t();
        else
            return f();
    }

    template<class Predicate, class Functor>
    [[nodiscard("[Warn]: Discarding a coroutine")]] auto co_if(Predicate&& pred, Functor&& func) noexcept(Functor()) {
        static_assert(std::predicate<Predicate>, "[Error]: Invalid predicate");
        using Handle = std::invoke_result<Functor>::type;
        if (pred())
            return func();
        else
            return Handle(nullptr);
    }
}
