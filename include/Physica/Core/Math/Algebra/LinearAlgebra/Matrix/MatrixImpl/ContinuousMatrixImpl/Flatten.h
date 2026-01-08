/*
 * Copyright 2023-2026 Weibo He.
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

#include "../ContinuousMatrix.h"

namespace Physica {
    template<Matrix T>
    class FlattenC<T> : public ContinuousVector<FlattenC<T>> {
        using This = FlattenC<T>;

        T& mat;
    public:
        using Base = ContinuousVector<This>;
    public:
        FlattenC(ContinuousMatrix<T>& mat_) : mat(mat_.getDerived()) {}
        FlattenC(const This&) = default;
        FlattenC(This&&) noexcept = default;
        ~FlattenC() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow() * mat.getCol(); }
        [[nodiscard]] auto data(this auto&& self) noexcept;
    };

    template<Matrix T>
    auto FlattenC<T>::data(this auto&& self) noexcept {
        return self.mat.data();
    }
}

namespace Physica {
    template<Matrix T>
    class Traits<FlattenC<T>> {
    public:
        using ScalarType = T::ScalarType;
        constexpr static size_t SizeAtCompile = T::RowAtCompile * T::ColAtCompile;

        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = true;
    };
}
