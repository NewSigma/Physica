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

#include "../CompactMatrix.cuh"

namespace Physica {
    template<class Derived>
    template<Matrix M>
    void device_obj<CompactMatrix<Derived>>::toHost(CompactMatrix<M>& obj) const {
        toHostAsync(obj);
        CUDAContext::getInstance().wait();
    }

    template<class Derived>
    template<Matrix M>
    void device_obj<CompactMatrix<Derived>>::toHostAsync(CompactMatrix<M>& obj) const {
        static_assert(std::same_as<ScalarType, typename M::ScalarType>, "[Error]: Incompatible ScalarType");
        obj.resize(Base::getRow(), Base::getCol());

        const size_t size = Base::getSize() * sizeof(T);
        if constexpr (M::getSizeAtCompile() != Dynamic)
            memcpy(obj.data(), data(), size);
        else if constexpr (Diffable<T>) {
            Base::getDerived().values().toHostAsync(obj.getDerived().values());
            Base::getDerived().grads().toHostAsync(obj.getDerived().grads());
        }
        else
            check(cudaMemcpyAsync(obj.data(), data(), size, cudaMemcpyKind::cudaMemcpyDeviceToHost, CUDAContext::getInstance()));
    }

    template<class Derived>
    __host__ __device__ auto device_obj<CompactMatrix<Derived>>::row(this auto&& self, size_t r) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        constexpr size_t ColAtCompile = Derived::getColAtCompile();
        constexpr bool IsMat1x1 = ColAtCompile == 1;
        if constexpr (IsMat1x1)
            return device_obj<CompactMatrixBlock<M, 1, 1>>(std::forward<Self>(self), r, 0);
        else {
            if constexpr (isRowMatrix)
                return device_obj<CompactMatrixBlock<M, 1, ColAtCompile>>(std::forward<Self>(self), r, 0, self.getCol());
            else
                return device_obj<LMatrixBlock<M, 1, ColAtCompile>>(std::forward<Self>(self), r, 0, self.getCol());
        }
    }

    template<class Derived>
    __host__ __device__ auto device_obj<CompactMatrix<Derived>>::col(this auto&& self, size_t c) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        constexpr size_t RowAtCompile = Derived::getRowAtCompile();
        constexpr bool IsMat1x1 = RowAtCompile == 1;
        if constexpr (IsMat1x1)
            return device_obj<CompactMatrixBlock<M, 1, 1>>(std::forward<Self>(self), 0, c);
        else {
            if constexpr (isColMatrix)
                return device_obj<CompactMatrixBlock<M, RowAtCompile, 1>>(std::forward<Self>(self), 0, self.getRow(), c);
            else
                return device_obj<LMatrixBlock<M, RowAtCompile, 1>>(std::forward<Self>(self), 0, self.getRow(), c);
        }
    }

    template<class Derived>
    template<size_t Row>
    __host__ __device__ auto device_obj<CompactMatrix<Derived>>::rows(this auto&& self, size_t fromRow, size_t rowCount) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<CompactMatrixBlock<M, Row, Derived::getColAtCompile()>>(std::forward<Self>(self), fromRow, rowCount, 0, self.getCol());
    }

    template<class Derived>
    template<size_t Row>
    __host__ __device__ auto device_obj<CompactMatrix<Derived>>::topRows(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<CompactMatrixBlock<M, Row, Derived::getColAtCompile()>>(std::forward<Self>(self), 0, to, 0, self.getCol());
    }

    template<class Derived>
    template<size_t Row>
    __host__ __device__ auto device_obj<CompactMatrix<Derived>>::bottomRows(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<CompactMatrixBlock<M, Row, Derived::getColAtCompile()>>(std::forward<Self>(self), from, self.getRow() - from, 0, self.getCol());
    }

    template<class Derived>
    template<size_t Col>
    __host__ __device__ auto device_obj<CompactMatrix<Derived>>::cols(this auto&& self, size_t fromCol, size_t colCount) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<CompactMatrixBlock<M, Derived::getRowAtCompile(), Col>>(std::forward<Self>(self), 0, self.getRow(), fromCol, colCount);
    }

    template<class Derived>
    template<size_t Col>
    __host__ __device__ auto device_obj<CompactMatrix<Derived>>::leftCols(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<CompactMatrixBlock<M, Derived::getRowAtCompile(), Col>>(std::forward<Self>(self), 0, self.getRow(), 0, to);
    }

    template<class Derived>
    template<size_t Col>
    __host__ __device__ auto device_obj<CompactMatrix<Derived>>::rightCols(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<CompactMatrixBlock<M, Derived::getRowAtCompile(), Col>>(std::forward<Self>(self), 0, self.getRow(), from, self.getCol() - from);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ auto device_obj<CompactMatrix<Derived>>::topLeftCorner(this auto&& self, size_t toRow, size_t toCol) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<CompactMatrixBlock<M, Row, Col>>(std::forward<Self>(self), 0, toRow, 0, toCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ auto device_obj<CompactMatrix<Derived>>::topLeftCorner(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<CompactMatrixBlock<M, Row, Col>>(std::forward<Self>(self), 0, to, 0, to);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ auto device_obj<CompactMatrix<Derived>>::topRightCorner(this auto&& self, size_t toRow, size_t fromCol) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<CompactMatrixBlock<M, Row, Col>>(std::forward<Self>(self), 0, toRow, fromCol, self.getRow() - fromCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ auto device_obj<CompactMatrix<Derived>>::bottomLeftCorner(this auto&& self, size_t fromRow, size_t toCol) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<CompactMatrixBlock<M, Row, Col>>(std::forward<Self>(self), fromRow, self.getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ auto device_obj<CompactMatrix<Derived>>::bottomRightCorner(this auto&& self, size_t fromRow, size_t fromCol) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<CompactMatrixBlock<M, Row, Col>>(std::forward<Self>(self), fromRow, self.getRow() - fromRow, fromCol, self.getCol() - fromCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ auto device_obj<CompactMatrix<Derived>>::bottomRightCorner(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<CompactMatrixBlock<M, Row, Col>>(std::forward<Self>(self), from, self.getRow() - from, from, self.getCol() - from);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ auto device_obj<CompactMatrix<Derived>>::block(this auto&& self, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<CompactMatrixBlock<M, Row, Col>>(std::forward<Self>(self), fromRow, rowCount, fromCol, colCount);
    }

    template<class Derived>
    void device_obj<CompactMatrix<Derived>>::zeros() {
        Base::getDerived().flatten().zeros();
    }

    template<class Derived>
    template<RNG R>
    void device_obj<CompactMatrix<Derived>>::random_uniform() {
        Base::getDerived().flatten().template random_uniform<R>();
    }

    template<class Derived>
    template<RNG R>
    void device_obj<CompactMatrix<Derived>>::random_normal() {
        Base::getDerived().flatten().template random_normal<R>();
    }

    template<class Derived>
    __host__ __device__ auto device_obj<CompactMatrix<Derived>>::data() noexcept {
        return Base::getDerived().data();
    }

    template<class Derived>
    __host__ __device__ auto device_obj<CompactMatrix<Derived>>::data() const noexcept {
        return Base::getDerived().data();
    }

    template<class Derived>
    __host__ __device__ auto device_obj<CompactMatrix<Derived>>::data_handle(this auto&& self) noexcept {
        return self.data();
    }

    template<class Derived>
    __host__ __device__ constexpr size_t device_obj<CompactMatrix<Derived>>::getRowStride() const noexcept {
        return Derived::isColMatrix() ? 1 : Base::getCol();
    }

    template<class Derived>
    __host__ __device__ constexpr size_t device_obj<CompactMatrix<Derived>>::getColStride() const noexcept {
        return Derived::isRowMatrix() ? 1 : Base::getRow();
    }

    template<class Derived>
    template<Matrix M>
    void CompactMatrix<Derived>::toDevice(device_obj<CompactMatrix<M>>& obj) const {
        toDeviceAsync(obj);
        if constexpr (!std::is_trivially_copy_constructible<M>::value)
            CUDAContext::getInstance().wait();
    }

    template<class Derived>
    template<Matrix M>
    void CompactMatrix<Derived>::toDeviceAsync(device_obj<CompactMatrix<M>>& obj) const {
        static_assert(std::is_same<T, typename M::ScalarType>::value,
                "[Error]: ScalarType inconsistent, additional buffer is necessary");
        obj.resize(Base::getRow(), Base::getCol());

        const size_t size = Base::getSize() * sizeof(T);
        if constexpr (M::getSizeAtCompile() != Dynamic)
            memcpy(obj.data(), data(), size);
        else if constexpr (Diffable<M>) {
            Base::getDerived().values().toDeviceAsync(obj.getDerived().values());
            Base::getDerived().grads().toDeviceAsync(obj.getDerived().grads());
        }
        else
            check(cudaMemcpyAsync(obj.data(), data(), size, cudaMemcpyKind::cudaMemcpyHostToDevice, CUDAContext::getInstance()));
    }
}
