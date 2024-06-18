/*
 * Copyright 2020-2023 WeiBo He.
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

namespace Physica::Core {
    template<ScalarOption Option>
    Scalar<Option> floor(const Scalar<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Scalar<Option>(::floorf(s.getTrivial()));
        else
            return Scalar<Option>(::floor(s.getTrivial()));
    }

    template<>
    Scalar<MultiPrecision> floor(const Scalar<MultiPrecision>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> ceil(const Scalar<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Scalar<Option>(::ceilf(s.getTrivial()));
        else
            return Scalar<Option>(::ceil(s.getTrivial()));
    }
}
