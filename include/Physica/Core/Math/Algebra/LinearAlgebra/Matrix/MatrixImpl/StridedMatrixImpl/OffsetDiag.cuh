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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/StridedVector.cuh"
#include "../StridedMatrix.cuh"

namespace Physica {
    template<Matrix M> requires(std::remove_cvref_t<M>::isStrided())
    class device_obj<OffsetDiag<M>> : public device_obj<StridedVector<OffsetDiag<M>>> {
        using host_obj = OffsetDiag<M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<StridedVector<host_obj>>;
        using Ref = add_device_obj<M>::type;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
        ssize_t offset;
    public:
        __host__ __device__ device_obj(Ref mat, ssize_t offset);
        device_obj(const This&) = default;
        device_obj(This&&) = default;
        ~device_obj() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t size) { assert(getLength() == size); }
        /* Getters */
        [[nodiscard]] __host__ __device__ auto&& getExpr(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle(this auto&& self) noexcept;
        [[nodiscard]] __host__ __device__ size_t getStride() const noexcept;
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept;
    };

    template<Matrix M> requires(std::remove_cvref_t<M>::isStrided())
    __host__ __device__ device_obj<OffsetDiag<M>>::device_obj(Ref mat, ssize_t offset) : mat(asStruct(mat)), offset(offset) {
        assert(getExpr().isSquare());
        assert(std::abs(offset) < getExpr().getRow());
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isStrided())
    __host__ __device__ auto&& device_obj<OffsetDiag<M>>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref>(self.mat.getDerived());
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isStrided())
    __host__ __device__ size_t device_obj<OffsetDiag<M>>::getLength() const noexcept {
        return getExpr().getRow() - std::abs(offset);
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isStrided())
    __host__ __device__ auto device_obj<OffsetDiag<M>>::data_handle(this auto&& self) noexcept {
        size_t r = self.offset < 0 ? -self.offset : 0;
        size_t c = self.offset > 0 ? self.offset : 0;
        return self.getExpr().data_handle() + r * self.getExpr().getRowStride() + c * self.getExpr().getColStride();
    }

    template<Matrix M> requires(std::remove_cvref_t<M>::isStrided())
    __host__ __device__ size_t device_obj<OffsetDiag<M>>::getStride() const noexcept {
        return getExpr().getRowStride() + getExpr().getColStride();
    }
}
