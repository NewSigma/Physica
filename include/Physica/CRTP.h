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
     * This class helps implementing CRTP(Curiously Recurring Template Pattern).
     * Typically, classes derived from \class CRTP are named with a 'Mixin' suffix", for example \class ScalarMixin, \class LayerMixin.
     *
     * Add host version since NVCC does not like SIMD
     */
    template<class T>
    class CRTP {
        using U = Traits<T>::Derived;
    public:
        [[nodiscard, gnu::always_inline, gnu::nodebug]] constexpr U& getDerived_host() noexcept { return *static_cast<U*>(this); }
        [[nodiscard, gnu::always_inline, gnu::nodebug]] constexpr const U& getDerived_host() const noexcept { return *static_cast<const U*>(this); }
        [[nodiscard, gnu::always_inline, gnu::nodebug]] constexpr U& getConstCastDerived_host() const noexcept { return *static_cast<U*>(const_cast<CRTP*>(this)); }

        [[nodiscard, gnu::always_inline, gnu::nodebug]] __host__ __device__ constexpr U& getDerived() noexcept { return *static_cast<U*>(this); }
        [[nodiscard, gnu::always_inline, gnu::nodebug]] __host__ __device__ constexpr const U& getDerived() const noexcept { return *static_cast<const U*>(this); }
        [[nodiscard, gnu::always_inline, gnu::nodebug]] __host__ __device__ constexpr U& getConstCastDerived() const noexcept { return *static_cast<U*>(const_cast<CRTP*>(this)); }
    protected:
        constexpr CRTP() = default;
        constexpr CRTP(const CRTP&) = default;
        constexpr CRTP(CRTP&&) noexcept = default;
        ~CRTP() = default;
        /* Operators */
        constexpr CRTP& operator=(const CRTP&) = default;
        constexpr CRTP& operator=(CRTP&&) noexcept = default;
    };
}
