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

#include <forward_list>
#include "CoDiff.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"

namespace Physica::Core {
    template<class Base>
    CoDiffNode<Base>::CoDiffNode(std::coroutine_handle<Promise> handle_) noexcept : Base(std::move(handle_.promise().obj)), handle(handle_) {}

    template<class Base>
    template<ReverseDiff T>
    CoDiffNode<Base>::CoDiffNode(T&& x) noexcept requires(!IsCoDiff<T>::value) {
        auto fn = [](T&& x) noexcept -> This {
            LazyDestroy<T&&> x_ = std::forward<T>(x);
            Base y;
            if constexpr (Scalar<T>) {
                auto result = co_yield Base(std::move(y));
                x_.reverse(result.grad());
            }
            else {
                static_assert(Vector<T> || Matrix<T>, "[Error]: Unexpected type T");
                y.resize(x_);
                x_.assign(y);
                auto result = co_yield std::move(y);
                x_.reverse(result.grads());
            }
        };
        fn(std::forward<T>(x)).swap(*this);
    }

    template<class Base>
    CoDiffNode<Base>::CoDiffNode(This&& other) noexcept : Base(std::move(other)), handle(other.handle) {
        other.handle = nullptr;
    }

    template<class Base>
    CoDiffNode<Base>::~CoDiffNode() {
        if (handle) {
            assert(!handle.done() && "[Error]: Unexpected resume, this is a bug");
            handle.promise().obj = Base(std::move(*this));
            handle.resume();
            assert(handle.done() && "[Error]: Invalid reverse diff");
            handle.destroy();
            handle = nullptr;
        }
    }

    template<class Base>
    void CoDiffNode<Base>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        std::swap(handle, obj.handle);
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
        else
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
        else
            return Handle(nullptr);
    }
}
