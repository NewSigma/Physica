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

#include "../StridedMatrix.cuh"
#include "Physica/Core/Exception/CUDA/CUDA.cuh"
#include "Physica/Core/Parallel/CUDAContext.cuh"

namespace Physica {
    template<class Derived>
    __host__ __device__ void device_obj<StridedMatrix<Derived>>::assign(this const auto& self, Matrix auto&& target) {
        using Self = std::remove_cvref_t<decltype(self)>;
        using Target = std::remove_cvref_t<decltype(target)>;
        target.assert_assign(self);
        if constexpr (host_obj::template isMinorCompact<Self, Target>() && IsHost()) {
            void* dst = target.data_handle();
            const void* src = self.data_handle();
            size_t dpitch = target.getMajorStride() * sizeof(T);
            size_t spitch = self.getMajorStride() * sizeof(T);
            size_t maxMajor = target.getMaxMajor();
            size_t maxMinor = target.getMaxMinor();
            assert(self.getMajorStride() >= maxMinor && target.getMajorStride() >= maxMinor && "[Error]: Rows must not overlap");
            check(cudaMemcpy2DAsync(dst, dpitch, src, spitch, maxMinor * sizeof(T), maxMajor, cudaMemcpyKind::cudaMemcpyDeviceToDevice, CUDAContext::getInstance()));
        }
        else
            self.device_obj<LValueMatrix<Derived>>::assign(target);
    }

    template<class Derived>
    __device__ void device_obj<StridedMatrix<Derived>>::assign(this const auto& self, Matrix auto&& target, instanceof_x<ThreadBlock> auto block) {
        self.device_obj<LValueMatrix<Derived>>::assign(std::forward<decltype(target)>(target), block);
    }

    template<class Derived>
    __host__ __device__ constexpr size_t device_obj<StridedMatrix<Derived>>::getRowStride() const noexcept {
        if constexpr (Derived::getRowStrideAtCompile() != Dynamic)
            return Derived::getRowStrideAtCompile();
        else
            return Base::getDerived().getRowStride();
    }

    template<class Derived>
    __host__ __device__ constexpr size_t device_obj<StridedMatrix<Derived>>::getColStride() const noexcept {
        if constexpr (Derived::getColStrideAtCompile() != Dynamic)
            return Derived::getColStrideAtCompile();
        else
            return Base::getDerived().getColStride();
    }

    template<class Derived>
    __host__ __device__ constexpr size_t device_obj<StridedMatrix<Derived>>::getMajorStride() const noexcept {
        if constexpr (Derived::isRowMatrix())
            return getRowStride();
        else
            return getColStride();
    }

    template<class Derived>
    __host__ __device__ constexpr size_t device_obj<StridedMatrix<Derived>>::getMinorStride() const noexcept {
        if constexpr (Derived::isRowMatrix())
            return getColStride();
        else
            return getRowStride();
    }

    template<class Derived>
    __host__ __device__ auto device_obj<StridedMatrix<Derived>>::data_handle() noexcept {
        return Base::getDerived().data_handle();
    }

    template<class Derived>
    __host__ __device__ auto device_obj<StridedMatrix<Derived>>::data_handle() const noexcept {
        return Base::getDerived().data_handle();
    }

    template<class Derived>
    __host__ __device__ auto device_obj<StridedMatrix<Derived>>::data_ptr(this auto&& self, size_t row, size_t col) noexcept {
        assert(row < self.getRow());
        assert(col < self.getCol());
        return self.data_handle() + row * self.getRowStride() + col * self.getColStride();
    }

    template<class Derived>
    __host__ __device__ auto device_obj<StridedMatrix<Derived>>::row(this auto&& self, size_t r) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<StridedMatrixBlock<M, 1, Derived::getColAtCompile()>>(std::forward<Self>(self), r, 1, 0, self.getCol());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<StridedMatrix<Derived>>::col(this auto&& self, size_t c) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<StridedMatrixBlock<M, Derived::getRowAtCompile(), 1>>(std::forward<Self>(self), 0, self.getRow(), c, 1);
    }

    template<class Derived>
    template<size_t Row>
    __host__ __device__ auto device_obj<StridedMatrix<Derived>>::rows(this auto&& self, size_t fromRow, size_t rowCount) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<StridedMatrixBlock<M, Row, Derived::getColAtCompile()>>(std::forward<Self>(self), fromRow, rowCount, 0, self.getCol());
    }

    template<class Derived>
    template<size_t Row>
    __host__ __device__ auto device_obj<StridedMatrix<Derived>>::topRows(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<StridedMatrixBlock<M, Row, Derived::getColAtCompile()>>(std::forward<Self>(self), 0, to, 0, self.getCol());
    }

    template<class Derived>
    template<size_t Row>
    __host__ __device__ auto device_obj<StridedMatrix<Derived>>::bottomRows(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<StridedMatrixBlock<M, Row, Derived::getColAtCompile()>>(std::forward<Self>(self), from, self.getRow() - from, 0, self.getCol());
    }

    template<class Derived>
    template<size_t Col>
    __host__ __device__ auto device_obj<StridedMatrix<Derived>>::cols(this auto&& self, size_t fromCol, size_t colCount) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<StridedMatrixBlock<M, Derived::getRowAtCompile(), Col>>(std::forward<Self>(self), 0, self.getRow(), fromCol, colCount);
    }

    template<class Derived>
    template<size_t Col>
    __host__ __device__ auto device_obj<StridedMatrix<Derived>>::leftCols(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<StridedMatrixBlock<M, Derived::getRowAtCompile(), Col>>(std::forward<Self>(self), 0, self.getRow(), 0, to);
    }

    template<class Derived>
    template<size_t Col>
    __host__ __device__ auto device_obj<StridedMatrix<Derived>>::rightCols(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<StridedMatrixBlock<M, Derived::getRowAtCompile(), Col>>(std::forward<Self>(self), 0, self.getRow(), from, self.getCol() - from);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ auto device_obj<StridedMatrix<Derived>>::topLeftCorner(this auto&& self, size_t toRow, size_t toCol) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<StridedMatrixBlock<M, Row, Col>>(std::forward<Self>(self), 0, toRow, 0, toCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ auto device_obj<StridedMatrix<Derived>>::topLeftCorner(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<StridedMatrixBlock<M, Row, Col>>(std::forward<Self>(self), 0, to, 0, to);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ auto device_obj<StridedMatrix<Derived>>::topRightCorner(this auto&& self, size_t toRow, size_t fromCol) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<StridedMatrixBlock<M, Row, Col>>(std::forward<Self>(self), 0, toRow, fromCol, self.getRow() - fromCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ auto device_obj<StridedMatrix<Derived>>::bottomLeftCorner(this auto&& self, size_t fromRow, size_t toCol) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<StridedMatrixBlock<M, Row, Col>>(std::forward<Self>(self), fromRow, self.getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ auto device_obj<StridedMatrix<Derived>>::bottomRightCorner(this auto&& self, size_t fromRow, size_t fromCol) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<StridedMatrixBlock<M, Row, Col>>(std::forward<Self>(self), fromRow, self.getRow() - fromRow, fromCol, self.getCol() - fromCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ auto device_obj<StridedMatrix<Derived>>::bottomRightCorner(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<StridedMatrixBlock<M, Row, Col>>(std::forward<Self>(self), from, self.getRow() - from, from, self.getCol() - from);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ auto device_obj<StridedMatrix<Derived>>::block(this auto&& self, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<StridedMatrixBlock<M, Row, Col>>(std::forward<Self>(self), fromRow, rowCount, fromCol, colCount);
    }

    template<class Derived>
    __host__ __device__ void device_obj<StridedMatrix<Derived>>::zeros() noexcept {
        if constexpr (host_obj::template isMinorCompact<Derived, Derived>() && IsHost()) {
            auto& self = Base::getDerived();
            const size_t maxMajor = self.getMaxMajor();
            const size_t maxMinor = self.getMaxMinor();
            assert(self.getMajorStride() >= maxMinor && "[Error]: Rows must not overlap");
            check(cudaMemset2DAsync(self.data_handle(), self.getMajorStride() * sizeof(T), 0, maxMinor * sizeof(T), maxMajor, CUDAContext::getInstance()));
        }
        else
            Base::zeros();
    }
}
