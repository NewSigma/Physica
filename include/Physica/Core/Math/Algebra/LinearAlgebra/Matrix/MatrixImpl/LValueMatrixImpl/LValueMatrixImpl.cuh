/*
 * Copyright 2024-2026 Weibo He.
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

#include "../LValueMatrix.cuh"

namespace Physica {
    template<class Derived>
    __host__ __device__ auto device_obj<LValueMatrix<Derived>>::operator=(Scalar auto x) -> device_obj<Derived>& {
        if (IsHost()) {
            auto func = [m_ = asStruct(Base::getDerived()), x] __device__() mutable {
                auto& m = m_.getDerived();
                const size_t maxMinor = m.getMaxMinor();
                const size_t major = blockIdx.y;
                const size_t minor = blockIdx.x * blockDim.x + threadIdx.x;
                if (minor < maxMinor)
                    m.refFromMajorMinor(major, minor) = x;
            };
            CUDAExecutor::launch<CUDADevAttr::DefaultThreadsPerBlock>(func, Base::makeKernelConfig());
        }
        else {
            for (size_t i = 0; i < Base::getMaxMajor(); ++i)
                for (size_t j = 0; j < Base::getMaxMinor(); ++j)
                    refFromMajorMinor(i, j) = x;
        }
        return Base::getDerived();
    }

    template<class Derived>
    __host__ __device__ void device_obj<LValueMatrix<Derived>>::operator+=(Scalar auto x) {
        auto& m = Base::getDerived();
        (m + x).assign(m);
    }

    template<class Derived>
    __host__ __device__ void device_obj<LValueMatrix<Derived>>::operator-=(Scalar auto x) {
        auto& m = Base::getDerived();
        (m - x).assign(m);
    }

    template<class Derived>
    __host__ __device__ void device_obj<LValueMatrix<Derived>>::operator*=(Scalar auto x) {
        auto& m = Base::getDerived();
        (m * x).assign(m);
    }

    template<class Derived>
    __host__ __device__ void device_obj<LValueMatrix<Derived>>::operator/=(Scalar auto x) {
        auto& m = Base::getDerived();
        (m / x).assign(m);
    }

    template<class Derived>
    __host__ __device__ device_obj<Derived>& device_obj<LValueMatrix<Derived>>::operator=(const Matrix auto& m) {
        auto& target = Base::getDerived();
        target.resize(m);
        m.assign(target);
        return target;
    }

    template<class Derived>
    __host__ __device__ void device_obj<LValueMatrix<Derived>>::operator+=(const Matrix auto& m) {
        m.assign_add(Base::getDerived());
    }

    template<class Derived>
    __host__ __device__ void device_obj<LValueMatrix<Derived>>::operator-=(const Matrix auto& m) {
        Base::getDerived() += -m;
    }

    template<class Derived>
    __device__ decltype(auto) device_obj<LValueMatrix<Derived>>::operator[](this auto&& self, size_t row, size_t col) {
        return *self.data_ptr(row, col);
    }

    template<class Derived>
    void device_obj<LValueMatrix<Derived>>::reverse(const Matrix auto& grad) const noexcept {
        using M = std::remove_cvref_t<decltype(grad)>;
        static_assert(std::same_as<typename ScalarType::GradType, typename M::ScalarType>, "[Error]: Inconsistent ScalarType");
        static_assert(isReverseDiff());
        assert(Base::getRow() == grad.getRow());
        assert(Base::getCol() == grad.getCol());
        grad.assign_add(Base::getConstCastDerived().grads());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<LValueMatrix<Derived>>::row(this auto&& self, size_t r) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<LMatrixBlock<M, 1, Dynamic>>(std::forward<Self>(self), r, 0, self.getCol());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<LValueMatrix<Derived>>::col(this auto&& self, size_t c) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<LMatrixBlock<M, Dynamic, 1>>(std::forward<Self>(self), 0, self.getRow(), c);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<LValueMatrix<Derived>>::rows(this auto&& self, size_t fromRow, size_t rowCount) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<RMatrixBlock<M>>(std::forward<Self>(self), fromRow, rowCount, 0, self.getCol());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<LValueMatrix<Derived>>::topRows(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<RMatrixBlock<M>>(std::forward<Self>(self), 0, to, 0, self.getCol());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<LValueMatrix<Derived>>::bottomRows(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<RMatrixBlock<M>>(std::forward<Self>(self), from, self.getRow() - from, 0, self.getCol());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<LValueMatrix<Derived>>::cols(this auto&& self, size_t fromCol, size_t colCount) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<RMatrixBlock<M>>(std::forward<Self>(self), 0, self.getRow(), fromCol, colCount);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<LValueMatrix<Derived>>::leftCols(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<RMatrixBlock<M>>(std::forward<Self>(self), 0, self.getRow(), 0, to);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<LValueMatrix<Derived>>::rightCols(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<RMatrixBlock<M>>(std::forward<Self>(self), 0, self.getRow(), from, self.getCol() - from);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<LValueMatrix<Derived>>::topLeftCorner(this auto&& self, size_t toRow, size_t toCol) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<RMatrixBlock<M>>(std::forward<Self>(self), 0, toRow, 0, toCol);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<LValueMatrix<Derived>>::topLeftCorner(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<RMatrixBlock<M>>(std::forward<Self>(self), 0, to, 0, to);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<LValueMatrix<Derived>>::topRightCorner(this auto&& self, size_t toRow, size_t fromCol) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<RMatrixBlock<M>>(std::forward<Self>(self), 0, toRow, fromCol, self.getRow() - fromCol);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<LValueMatrix<Derived>>::bottomLeftCorner(this auto&& self, size_t fromRow, size_t toCol) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<RMatrixBlock<M>>(std::forward<Self>(self), fromRow, self.getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<LValueMatrix<Derived>>::bottomRightCorner(this auto&& self, size_t fromRow, size_t fromCol) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<RMatrixBlock<M>>(std::forward<Self>(self), fromRow, self.getRow() - fromRow, fromCol, self.getCol() - fromCol);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<LValueMatrix<Derived>>::bottomRightCorner(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<RMatrixBlock<M>>(std::forward<Self>(self), from, self.getRow() - from, from, self.getCol() - from);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<LValueMatrix<Derived>>::block(this auto&& self, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<RMatrixBlock<M>>(std::forward<Self>(self), fromRow, rowCount, fromCol, colCount);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<LValueMatrix<Derived>>::diag(this auto&& self) noexcept {
        assert(self.isSquare());
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<DiagVectorL<M>>(std::forward<Self>(self));
    }

    template<class Derived>
    __host__ __device__ auto device_obj<LValueMatrix<Derived>>::flatten() {
        return device_obj<FlattenL<Derived>>(Base::getDerived());
    }

    template<class Derived>
    __host__ __device__ const auto device_obj<LValueMatrix<Derived>>::flatten() const {
        return device_obj<FlattenL<Derived>>(Base::getDerived());
    }

    template<class Derived>
    void device_obj<LValueMatrix<Derived>>::zero_grad() noexcept {
        Base::getDerived().grads().zeros();
    }

    template<class Derived>
    template<RNG R>
    void device_obj<LValueMatrix<Derived>>::random_uniform() {
        Derived::template random_uniform<R>(Base::getRow(), Base::getCol()).toDeviceAsync(Base::getDerived());
    }

    template<class Derived>
    template<RNG R>
    void device_obj<LValueMatrix<Derived>>::random_normal() {
        Derived::template random_normal<R>(Base::getRow(), Base::getCol()).toDeviceAsync(Base::getDerived());
    }

    template<class Derived>
    template<RNG R>
    void device_obj<LValueMatrix<Derived>>::random_any(auto& distribution) {
        Derived::template random_any<R>(Base::getRow(), Base::getCol(), distribution).toDeviceAsync(Base::getDerived());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<LValueMatrix<Derived>>::data_ptr(this auto&& self, size_t row, size_t col) noexcept {
        assert(row < self.getRow());
        assert(col < self.getCol());
        return self.getDerived().data_ptr(row, col);
    }

    template<class Derived>
    __device__ decltype(auto) device_obj<LValueMatrix<Derived>>::refFromMajorMinor(this auto&& self, size_t major, size_t minor) noexcept {
        assert(major < self.getMaxMajor());
        assert(minor < self.getMaxMinor());
        const size_t r = MatrixMajor::rowFromMajorMinor<Derived>(major, minor);
        const size_t c = MatrixMajor::colFromMajorMinor<Derived>(major, minor);
        return self[r, c];
    }
}
