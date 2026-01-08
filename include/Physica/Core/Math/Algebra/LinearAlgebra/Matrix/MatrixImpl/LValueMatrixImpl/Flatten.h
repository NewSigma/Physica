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

#include "../LValueMatrix.h"

namespace Physica {
    template<Matrix T>
    class FlattenL<T> : public LValueVector<FlattenL<T>> {
        using This = FlattenL<T>;
        using Base = LValueVector<This>;

        T& mat;
    public:
        FlattenL(LValueMatrix<T>& mat_) : mat(mat_.getDerived()) {}
        FlattenL(const This&) = default;
        FlattenL(This&&) noexcept = default;
        ~FlattenL() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow() * mat.getCol(); }
        [[nodiscard]] auto data_ptr(this auto&&, size_t index) noexcept;
    };

    template<Matrix T>
    auto FlattenL<T>::data_ptr(this auto&& self, size_t index) noexcept {
        const size_t major = index / self.mat.getMaxMinor();
        const size_t minor = index % self.mat.getMaxMinor();
        const size_t row = MatrixOption::rowFromMajorMinor<T>(major, minor);
        const size_t col = MatrixOption::colFromMajorMinor<T>(major, minor);
        return self.mat.data_ptr(row, col);
    }
}

namespace Physica {
    template<Matrix T>
    class Traits<FlattenL<T>> {
    public:
        using ScalarType = T::ScalarType;
        constexpr static size_t SizeAtCompile = T::RowAtCompile * T::ColAtCompile;

        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = false;
    };
}
