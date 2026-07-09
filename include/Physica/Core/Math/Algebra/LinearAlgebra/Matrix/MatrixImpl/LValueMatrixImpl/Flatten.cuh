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

#include "Flatten.h"

namespace Physica {
    template<Matrix M> requires (std::remove_cvref_t<M>::isLValueMatrix() && !std::remove_cvref_t<M>::isCompact())
    class device_obj<Flatten<M>> : public device_obj<LValueVector<Flatten<M>>> {
        using host_obj = Flatten<M>;
        using This = device_obj<host_obj>;

        const device_obj<M>& mat;
    public:
        using Base = device_obj<LValueVector<host_obj>>;
    public:
        __host__ __device__ device_obj(const device_obj<LValueMatrix<M>>& mat_) : mat(mat_) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        using Base::resize;
        __host__ __device__ void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return mat.getRow() * mat.getCol(); }
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&& self, size_t index) noexcept;
    };

    template<Matrix M> requires (std::remove_cvref_t<M>::isLValueMatrix() && !std::remove_cvref_t<M>::isCompact())
    __host__ __device__ auto device_obj<Flatten<M>>::data_ptr(this auto&& self, size_t index) noexcept {
        auto&& mat = std::forward<decltype(self)>(self).mat;
        const size_t major = index / mat.getMaxMinor();
        const size_t minor = index % mat.getMaxMinor();
        const size_t row = MatrixMajor::rowFromMajorMinor<M>(major, minor);
        const size_t col = MatrixMajor::colFromMajorMinor<M>(major, minor);
        return mat.data_ptr(row, col);
    }
}
