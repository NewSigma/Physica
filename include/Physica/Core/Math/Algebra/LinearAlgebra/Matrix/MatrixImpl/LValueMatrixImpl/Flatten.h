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
    template<Matrix M>
    class FlattenL<M> : public LValueVector<FlattenL<M>> {
        using This = FlattenL<M>;
        using Base = LValueVector<This>;

        LazyDestroy<M> mat;
    public:
        FlattenL(M&& mat_) : mat(std::forward<M>(mat_)) {}
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

        [[nodiscard]] auto values(this auto&&) noexcept;
        [[nodiscard]] auto grads(this auto&& self) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow() * mat.getCol(); }
        [[nodiscard]] auto data_ptr(this auto&&, size_t index) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept;
    };

    template<Matrix M>
    auto FlattenL<M>::values(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat).values().flatten();
    }

    template<Matrix M>
    auto FlattenL<M>::grads(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat).grads().flatten();
    }

    template<Matrix M>
    auto FlattenL<M>::data_ptr(this auto&& self, size_t index) noexcept {
        const size_t major = index / self.mat.getMaxMinor();
        const size_t minor = index % self.mat.getMaxMinor();
        const size_t row = MatrixMajor::rowFromMajorMinor<M>(major, minor);
        const size_t col = MatrixMajor::colFromMajorMinor<M>(major, minor);
        return self.mat.data_ptr(row, col);
    }

    template<Matrix M>
    __host__ __device__ consteval size_t FlattenL<M>::getSizeAtCompile() noexcept {
        using M1 = std::remove_reference_t<M>;
        return M1::RowAtCompile * M1::ColAtCompile;
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<FlattenL<M>> {
    public:
        using ScalarType = std::remove_reference_t<M>::ScalarType;
        constexpr static bool FastAssign = false;
    };
}
