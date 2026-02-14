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
    template<class Base>
    class DiffCoro : public Base {
        static_assert(ReverseDiff<Base>, "[Error]: DiffCoro save compute graph for reverse diffable objects");
        static_assert(!is_codiff<Base>::value, "[Error]: Nested DiffCoro is not allowed");
        static_assert(std::is_object<Base>::value, "[Error]: Must save the return by value");
        using This = DiffCoro<Base>;

        class Promise {
            struct suspend_yield : public std::suspend_always {
                Promise* p;

                explicit suspend_yield(Promise* p) : p(p) {}

                static void await_suspend(std::coroutine_handle<>) noexcept {
                    // In LLVM, await_suspend is implemented as an intrinsic. Making it static would enable more optimizations.
                }

                [[nodiscard]] Base& await_resume() const noexcept { return p->obj; }
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
            auto get_return_object() noexcept { return std::coroutine_handle<Promise>::from_promise(*this); };
            static auto get_return_object_on_allocation_failure() noexcept { return nullptr; }
            std::suspend_never initial_suspend() noexcept { return {}; }
            void await_transform(auto&&) noexcept = delete("[Error]: Differential coroutine must yield values");
            std::suspend_never final_suspend() noexcept { return {}; }
            suspend_yield yield_value(auto&& arg) noexcept;
            void return_void() noexcept {}
            [[noreturn]] void unhandled_exception() { throw; }

            void swap(This& __restrict obj_) noexcept { obj.swap(obj_); }
        };
    public:
        using promise_type = Promise;
        using typename Base::ScalarType;
    private:
        std::coroutine_handle<Promise> handle = nullptr;
    public:
        DiffCoro() = default;
        DiffCoro(std::nullptr_t) noexcept {}
        DiffCoro(std::coroutine_handle<Promise> handle_) noexcept;
        template<ReverseDiff T>
        DiffCoro(T&& x) noexcept requires(!is_codiff<T>::value); // Codiff should delegate to the move constructor
        DiffCoro(const This& other) = delete;
        DiffCoro(This&& other) noexcept;
        ~DiffCoro();
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&& obj) noexcept;
        using Base::operator=;
        /* Operations */
        void reverse_final() noexcept;
        void reverse_final(auto&& x) noexcept;
        void swap(This& __restrict obj) noexcept;
    private:
        void reverse_impl() noexcept;
        [[nodiscard]] static This compute(ReverseDiff auto&& expr) noexcept;
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
    template<class Base>
    struct formatter<Physica::DiffCoro<Base>, char> {
        constexpr auto parse(std::format_parse_context& ctx) {
            return ctx.begin();
        }

        static auto format(const Physica::DiffCoro<Base>& obj, auto& ctx) {
            return formatter<Base, char>{}.format(obj, ctx);
        }
    };
}

#include "DiffCoroImpl.h"
