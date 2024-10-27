/*
 * Copyright 2021-2024 Weibo He.
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

#include <cstddef>
#include "Macro.h"

namespace Physica {
    /**
     * This class helps implementing CRTP.
     */
    template<class T>
    class CRTPBase {
        using This = CRTPBase<T>;
        using U = typename Traits<T>::Derived;
    public:
        [[nodiscard]] __host__ __device__ U& getDerived() noexcept { return *static_cast<U*>(this); }
        [[nodiscard]] __host__ __device__ const U& getDerived() const noexcept { return *static_cast<const U*>(this); }
        [[nodiscard]] __host__ __device__ U& getConstCastDerived() const noexcept { return *static_cast<U*>(const_cast<CRTPBase*>(this)); }
    protected:
        CRTPBase() = default;
        CRTPBase(const This&) = default;
        CRTPBase(This&&) noexcept = default;
        ~CRTPBase() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
    };
}
