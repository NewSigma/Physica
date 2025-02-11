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

#include <cstdlib>
#include <exception>
#include <coroutine>
#include "CRTPBase.h"

namespace Physica {
    /**
     * \class CRCoro is Curiously Recurring Coroutine
     *
     * Allow a normal function become a trivial coroutine without cost
     */
    template<class T>
    class CRCoro : public CRTPBase<CRCoro<T>> {
        static_assert(std::is_object<T>::value, "[Error]: Must save the return by value");
        using This = CRCoro<T>;
        using Base = CRTPBase<CRCoro<T>>;

        struct RValueWrapper {
            T* p;

            ~RValueWrapper() {
                std::coroutine_handle<T>::from_promise(*p).destroy();
            }

            operator T&&() const noexcept { return std::move(*p); }
        };
    public:
        using promise_type = T;
    public:
        auto get_return_object() noexcept { return RValueWrapper(&Base::getDerived()); }
        std::suspend_never initial_suspend() noexcept { return {}; }
        std::suspend_always final_suspend() noexcept { return {}; }
        void return_value(T&& x) noexcept { Base::getDerived() = std::move(x); }
        void unhandled_exception() { throw std::current_exception(); }
    protected:
        CRCoro() = default;
        CRCoro(const This&) = default;
        CRCoro(This&&) noexcept = default;
        ~CRCoro() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
    };
}

namespace Physica {
    template<class T>
    class Traits<CRCoro<T>> {
    public:
        using Derived = T;
    };
}
