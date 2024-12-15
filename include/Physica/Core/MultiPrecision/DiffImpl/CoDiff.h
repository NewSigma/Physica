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

#include <cassert>
#include <coroutine>
#include <exception>
#include <forward_list>
#include "Physica/Core/MultiPrecision/Scalar.h"

namespace Physica::Core {
    template<class T>
    struct IsCoDiffNode {
        constexpr static bool value = false;
    };

    template<class T>
    struct IsCoDiffNode<CoDiffNode<T>> {
        constexpr static bool value = false;
    };

    template<class T> requires(std::is_reference<T>::value)
    using LazyReverse = std::conditional<std::is_rvalue_reference<T>::value, std::remove_reference_t<T>, T&>::type;

    template<>
    class CoDiffNode<void> {
        using This = CoDiffNode<void>;
        class Impl {
        public:
            auto get_return_object() { return std::coroutine_handle<Impl>::from_promise(*this); };
            std::suspend_never initial_suspend() noexcept { return {}; }
            std::suspend_always final_suspend() noexcept { return {}; }
            void return_void() noexcept {}
            void unhandled_exception() { throw std::current_exception(); }
        };
    public:
        using promise_type = Impl;
    private:
        std::coroutine_handle<Impl> handle = nullptr;
    public:
        CoDiffNode() = default;
        CoDiffNode(std::coroutine_handle<Impl> handle_) noexcept : handle(handle_) {}
        CoDiffNode(const This&) = delete;
        CoDiffNode(This&& other) noexcept : handle(other.handle) { other.handle = nullptr; }
        ~CoDiffNode() {
            if (handle) {
                assert(!handle.done() && "[Error]: Unexpected resume, this is a bug");
                handle.resume();
                assert(handle.done() && "[Error]: Invalid reverse diff");
                handle.destroy();
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
    };

    template<class Base>
    class CoDiffNode : public Base {
        static_assert(ReverseDiff<Base>, "[Error]: CoDiffNode save compute graph for reverse diffable objects");
        static_assert(!IsCoDiffNode<Base>::value, "[Error]: Nested CoDiffNode is not allowed");
        static_assert(std::is_object<Base>::value, "[Error]: Must save the return by value");
        using This = CoDiffNode<Base>;

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
        CoDiffNode() = default;
        CoDiffNode(std::coroutine_handle<Impl> handle_) noexcept;
        template<ReverseDiff T>
        CoDiffNode(T&& x) noexcept;
        CoDiffNode(const This&) = delete;
        CoDiffNode(This&& other) noexcept;
        ~CoDiffNode();
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(This& __restrict obj) noexcept;
    };

    template<class Predicate, class Operation, class Functor>
    [[nodiscard("[Warn]: Discarding a coroutine")]] auto co_for(Predicate&& pred, Operation&& op, Functor&& func) noexcept(noexcept(func()));

    template<class Predicate, class FuncTrue, class FuncFalse>
    [[nodiscard("[Warn]: Discarding a coroutine")]] auto co_if(Predicate&& pred, FuncTrue&& t, FuncFalse&& f) noexcept(noexcept(t()) && noexcept(f()));

    template<class Predicate, class Functor>
    [[nodiscard("[Warn]: Discarding a coroutine")]] auto co_if(Predicate&& pred, Functor&& func) noexcept(noexcept(func()));
}

#include "CoDiffImpl.h"
