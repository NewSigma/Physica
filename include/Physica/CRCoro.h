/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <cstdlib>
#include <utility>
#include "CRTP.h"
#include "Physica/Core/Utils/Suspend.h"

namespace Physica {
    /**
     * \class CRCoro is Curiously Recurring Coroutine
     *
     * Allow a normal function become a trivial coroutine without cost
     */
    template<class T>
    class CRCoro : public CRTP<CRCoro<T>> {
        static_assert(std::is_object<T>::value, "[Error]: Must save the return by value");
        using This = CRCoro<T>;
        using Base = CRTP<This>;
    public:
        using promise_type = T;
    public:
        ~CRCoro() = default;
        /* Operations */
        auto get_return_object() noexcept {
            struct RValueWrapper {
                This* p;

                operator T&&() const noexcept { return std::move(p->getDerived()); }
            };
            return RValueWrapper(this);
        }
        auto initial_suspend() noexcept { return suspend_never{}; }
        void await_transform(auto&&) noexcept = delete("[Error]: CRCoro never suspends");
        auto final_suspend() noexcept { return suspend_never{}; }
        void return_value(auto&& x) noexcept { Base::getDerived() = std::forward<decltype(x)>(x); }
        [[noreturn]] void unhandled_exception() { throw; }
    protected:
        CRCoro() = default;
        CRCoro(const This&) = default;
        CRCoro(This&&) noexcept = default;
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
