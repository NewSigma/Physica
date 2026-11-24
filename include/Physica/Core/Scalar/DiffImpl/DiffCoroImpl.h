/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

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
        assert(handle != nullptr);
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
        if (handle)
            handle.promise().listen(*this);
        if (obj.handle)
            obj.handle.promise().listen(obj);
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
}
