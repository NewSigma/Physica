/*
 * Copyright 2022-2025 Weibo He.
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

#include "../RValueMatrix.h"

namespace Physica {
    template<Matrix T>
    class FlattenR<T> : public RValueVector<FlattenR<T>> {
        using This = FlattenR<T>;

        const T& mat;
    public:
        using Base = RValueVector<FlattenR<T>>;
        using typename Base::ScalarType;
    public:
        FlattenR(const T& mat_) : mat(mat_) {}
        FlattenR(const This&) = default;
        FlattenR(This&&) noexcept = default;
        ~FlattenR() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t index) const;
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow() * mat.getCol(); }
    };

    template<Matrix T>
    auto FlattenR<T>::calc(size_t index) const -> ScalarType {
        const size_t major = index / mat.getMaxMinor();
        const size_t minor = index % mat.getMaxMinor();
        return mat.calcFromMajorMinor(major, minor);
    }
}

namespace Physica {
    template<Matrix T>
    class Traits<FlattenR<T>> {
    public:
        using ScalarType = T::ScalarType;
        constexpr static size_t SizeAtCompile = T::RowAtCompile * T::ColAtCompile;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = false;
    };
}
