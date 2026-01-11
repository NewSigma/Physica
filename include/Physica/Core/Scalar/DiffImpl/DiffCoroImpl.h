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

#include <forward_list>
#include "DiffCoro.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"

namespace Physica {
    template<class Base>
    auto DiffCoro<Base>::Promise::yield_value(auto&& arg) noexcept -> suspend_yield {
        if constexpr (Scalar<decltype(arg)>)
            obj = Base(arg.value());
        else
            obj = Base(std::forward<decltype(arg)>(arg).values());
        return suspend_yield{};
    }

    template<class Base>
    DiffCoro<Base>::DiffCoro(std::coroutine_handle<Promise> handle_) noexcept : Base(std::move(handle_.promise().obj)), handle(handle_) {
        /**
         * Reference:
         * [1] GH123347; https://github.com/llvm/llvm-project/issues/123347
         */
    #if defined(__clang__) && (__clang_major__ <= 20)
        asm volatile("" : : "r,m"(handle) : "memory");
    #endif
    }
    /**
     * Lazily compute expression \tparam T and construct compute graph from the result
     */
    template<class Base>
    template<ReverseDiff T>
    DiffCoro<Base>::DiffCoro(T&& x) noexcept requires(!is_codiff<T>::value) : DiffCoro(compute(std::forward<T>(x))) {}

    template<class Base>
    DiffCoro<Base>::DiffCoro(This&& other) noexcept : Base(static_cast<Base&&>(other)), handle(other.handle) {
        other.handle = nullptr;
    }

    template<class Base>
    DiffCoro<Base>::~DiffCoro() {
        reverse_impl();
    }

    template<class Base>
    auto DiffCoro<Base>::operator=(This&& obj) noexcept -> This& {
        swap(obj);
        return *this;
    }

    template<class Base>
    void DiffCoro<Base>::reverse_final() noexcept {
        assert(handle != nullptr && "[Error]: Reverse has been finalized");
        Base::reverse();
        reverse_impl();
    }

    template<class Base>
    void DiffCoro<Base>::reverse_final(auto&& x) noexcept {
        assert(handle != nullptr && "[Error]: Reverse has been finalized");
        Base::reverse(std::forward<decltype(x)>(x));
        reverse_impl();
    }

    template<class Base>
    void DiffCoro<Base>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        std::swap(handle, obj.handle);
    }

    template<class Base>
    void DiffCoro<Base>::reverse_impl() noexcept {
        if (handle) {
            assert(!handle.done() && "[Error]: Unexpected resume, this is a bug");
            if constexpr (Scalar<Base>) {
                if (Base::grad().isZero()) {
                    handle.destroy();
                    handle = nullptr;
                    return;
                }
            }
            handle.promise().obj = Base(static_cast<Base&&>(*this));
            handle.resume();
            handle = nullptr;
        }
    }

    template<class Base>
    template<ReverseDiff T>
    auto DiffCoro<Base>::compute(T&& expr) noexcept -> This {
        static_assert(!std::same_as<std::remove_cvref_t<T>, Base>, "[Error]: Not a expression");
        static_assert(Vector<T> || Matrix<T>, "[Error]: Not a expression");

        LazyDestroy<T&&> expr_ = std::forward<T>(expr);
        auto result = co_yield expr_.values();
        expr_.reverse(result.values(), result.grads());
    }

    template<class Predicate, class Operation, class Functor>
    [[nodiscard("[Warn]: Discarding a coroutine")]] auto co_for(Predicate&& pred, Operation&& op, Functor&& func) noexcept(noexcept(func())) {
        using Handle = std::invoke_result<Functor>::type;
        static_assert(std::is_class<Handle>::value, "[Error]: Invalid handle type");
        static_assert(!std::is_trivially_destructible<Handle>::value, "[Error]: Invalid handle type");
        static_assert(std::predicate<Predicate>, "[Error]: Invalid predicate");

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
        static_assert(std::is_class<Handle1>::value, "[Error]: Invalid handle type");
        static_assert(!std::is_trivially_destructible<Handle1>::value, "[Error]: Invalid handle type");
        static_assert(std::predicate<Predicate>, "[Error]: Invalid predicate");

        if (pred())
            return t();
        return f();
    }

    template<class Predicate, class Functor>
    [[nodiscard("[Warn]: Discarding a coroutine")]] auto co_if(Predicate&& pred, Functor&& func) noexcept(noexcept(func())) {
        using Handle = std::invoke_result<Functor>::type;
        static_assert(std::is_class<Handle>::value, "[Error]: Invalid handle type");
        static_assert(!std::is_trivially_destructible<Handle>::value, "[Error]: Invalid handle type");
        static_assert(std::predicate<Predicate>, "[Error]: Invalid predicate");

        if (pred())
            return func();
        return Handle(nullptr);
    }
}
