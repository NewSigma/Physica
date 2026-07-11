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
    /**
     * Lazily compute expression and construct compute graph from the result
     */
    template<class T>
    DiffCoro<T>::DiffCoro(ReverseDiff auto&& expr) noexcept requires(!is_codiff<decltype(expr)>::value)
            : DiffCoro(compute(std::forward<decltype(expr)>(expr))) {}

    template<class T>
    DiffCoro<T>::DiffCoro(Promise& p) noexcept : handle(std::coroutine_handle<Promise>::from_promise(p)) {
        p.listen(*this);
    }

    template<class T>
    DiffCoro<T>::DiffCoro(This&& other) noexcept
            : Base(static_cast<T&&>(other))
            , handle(std::exchange(other.handle, nullptr)) {
        handle.promise().listen(*this);
    }

    template<class T>
    DiffCoro<T>::~DiffCoro() {
        reverse_impl();
    }

    template<class T>
    auto DiffCoro<T>::operator=(This&& obj) noexcept -> This& {
        swap(obj);
        return *this;
    }

    template<class T>
    void DiffCoro<T>::reverse_final(auto&&... args) noexcept {
        assert(handle != nullptr && "[Error]: Reverse has been finalized");
        Base::reverse(std::forward<decltype(args)>(args)...);
        reverse_impl();
    }

    template<class T>
    void DiffCoro<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        std::swap(handle, obj.handle);
    }

    template<class T>
    void DiffCoro<T>::reverse_impl() noexcept {
        if (handle) {
            assert(!handle.done() && "[Error]: Unexpected resume, this is a bug");
            handle.resume();
            handle = nullptr;
        }
    }

    template<class T>
    auto DiffCoro<T>::compute(ReverseDiff auto&& expr) noexcept -> This {
        using Expr = decltype(expr);
        static_assert(!std::same_as<T, std::remove_cvref_t<Expr>>, "[Error]: Not a expression");
        static_assert(Vector<Expr> || Matrix<Expr>, "[Error]: Not a expression");

        decltype(auto) expr_ = decay_rvalue(std::forward<Expr>(expr));
        auto& result = co_yield expr_.values();
        expr_.reverse(result.values(), result.grads());
    }

    template<class T>
    void DiffCoro<T>::Promise::listen(DiffCoro<T>& node) noexcept {
        pObj = &node; // Other nodes may take ownership of the promise
    }

    template<class T>
    auto DiffCoro<T>::Promise::get_return_object() noexcept {
        return DiffCoro<T>(*this);
    };

    template<class T>
    auto& DiffCoro<T>::Promise::yield_value(auto&& arg) noexcept {
        using Arg = decltype(arg);
        if constexpr (Scalar<Arg> || Packet<Arg>)
            *pObj = T(arg.value());
        else
            *pObj = T(std::forward<Arg>(arg).values());
        return *this;
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
