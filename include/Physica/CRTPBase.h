/*
 * Copyright 2021-2026 Weibo He.
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

#include "Macro.h"

namespace Physica {
    /**
     * This class helps implementing CRTP.
     *
     * Add host version since NVCC does not like SIMD
     */
    template<class T>
    class CRTPBase {
        using This = CRTPBase<T>;
        using U = Traits<T>::Derived;
    public:
        [[nodiscard, gnu::always_inline, gnu::nodebug]] constexpr U& getDerived_host() noexcept { return *static_cast<U*>(this); }
        [[nodiscard, gnu::always_inline, gnu::nodebug]] constexpr const U& getDerived_host() const noexcept { return *static_cast<const U*>(this); }
        [[nodiscard, gnu::always_inline, gnu::nodebug]] constexpr U& getConstCastDerived_host() const noexcept { return *static_cast<U*>(const_cast<CRTPBase*>(this)); }

        [[nodiscard, gnu::always_inline, gnu::nodebug]] __host__ __device__ constexpr U& getDerived() noexcept { return *static_cast<U*>(this); }
        [[nodiscard, gnu::always_inline, gnu::nodebug]] __host__ __device__ constexpr const U& getDerived() const noexcept { return *static_cast<const U*>(this); }
        [[nodiscard, gnu::always_inline, gnu::nodebug]] __host__ __device__ constexpr U& getConstCastDerived() const noexcept { return *static_cast<U*>(const_cast<CRTPBase*>(this)); }
    protected:
        constexpr CRTPBase() = default;
        constexpr CRTPBase(const This&) = default;
        constexpr CRTPBase(This&&) noexcept = default;
        ~CRTPBase() = default;
        /* Operators */
        constexpr This& operator=(const This&) = default;
        constexpr This& operator=(This&&) noexcept = default;
    };
}
