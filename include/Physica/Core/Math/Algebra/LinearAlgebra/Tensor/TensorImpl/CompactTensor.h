/*
 * Copyright 2026 Weibo He.
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

#include "StridedTensor.h"

namespace Physica {
    template<class Derived>
    class CompactTensor : public StridedTensor<Derived> {
        using Base = StridedTensor<Derived>;
        using This = CompactTensor<Derived>;
    public:
        using typename Base::IndexType;
    public:
        ~CompactTensor() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        using Base::operator=;
        /* Operations */
        [[nodiscard]] auto fiber(this auto&&, IndexVar auto...) noexcept;
        /* Getters */
        [[nodiscard]] auto data(this auto&& self) noexcept;
        [[nodiscard]] auto data_handle() noexcept;
        [[nodiscard]] auto data_handle() const noexcept;
        [[nodiscard]] constexpr auto getStrides() const noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isCompact() noexcept { return true; }
    protected:
        CompactTensor() = default;
        CompactTensor(const This&) = default;
        CompactTensor(This&&) noexcept = default;
    };
}

namespace Physica {
    template<class Derived>
    class Traits<CompactTensor<Derived>> : public Traits<Derived> {};
}

#include "CompactTensorImpl/CompactTensorImpl.h"
