/*
 * Copyright 2022-2026 Weibo He.
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

#include "../RValueMatrix.cuh"
#include "Physica/PlainStruct.h"
#include "Flatten.cuh" // IWYU pragma: export

namespace Physica {
    template<class Derived>
    template<Matrix M>
    __host__ __device__ void device_obj<RValueMatrix<Derived>>::assign(M&& target) const {
        target.assert_assign(Base::getDerived());
        if (IsHost()) {
            auto func = [source_ = asStruct(Base::getDerived()), target_ = asStruct(target)] __device__() mutable {
                const auto& source = source_.getDerived();
                auto& target = target_.getDerived();
                const size_t maxMinor = target.getMaxMinor();
                const size_t minor = blockIdx.x * blockDim.x + threadIdx.x;
                if (minor < maxMinor) {
                    const size_t major = blockIdx.y;
                    const size_t r = target.rowFromMajorMinor(major, minor);
                    const size_t c = target.colFromMajorMinor(major, minor);
                    target[r, c] = source.calc(r, c);
                }
            };
            CUDAExecutor::launch<MaxThreadsPerBlock>(func, target.makeKernelConfig());
        }
        else if constexpr (IsDevice()) {
            const size_t maxMajor = target.getMaxMajor();
            const size_t maxMinor = target.getMaxMinor();
            for (size_t major = 0; major < maxMajor; ++major) {
                for (size_t minor = 0; minor < maxMinor; ++minor) {
                    const size_t r = target.rowFromMajorMinor(major, minor);
                    const size_t c = target.colFromMajorMinor(major, minor);
                    target[r, c] = Base::getDerived().calc(r, c);
                }
            }
        }
    }

    template<class Derived>
    __host__ __device__ void device_obj<RValueMatrix<Derived>>::assign_add(Matrix auto&& target) const {
        const auto& x = Base::getDerived();
        target.assert_assign(x);
        target = target + x;
    }

    template<class Derived>
    __host__ __device__ void device_obj<RValueMatrix<Derived>>::assert_assign(const Matrix auto& source) const noexcept {
        static_assert_assign(source);
        assert(getRow() == source.getRow() && "[Error]: Dimensions do not match");
        assert(getCol() == source.getCol() && "[Error]: Dimensions do not match");
    }

    template<class Derived>
    __host__ __device__ void device_obj<RValueMatrix<Derived>>::resize(const Matrix auto& m, auto&&... args) {
        resize(m.getRow(), m.getCol(), std::forward<decltype(args)>(args)...);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::resize(size_t r, size_t c, auto&&... args) {
        return Base::getDerived().resize(r, c, std::forward<decltype(args)>(args)...);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::row(size_t r) noexcept {
        return RowVector(Base::getDerived(), r, 0, getCol());
    }

    template<class Derived>
    __host__ __device__ const auto device_obj<RValueMatrix<Derived>>::row(size_t r) const noexcept {
        return RowVector(Base::getConstCastDerived(), r, 0, getCol());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::col(size_t c) noexcept {
        return ColVector(Base::getDerived(), 0, getRow(), c);
    }

    template<class Derived>
    __host__ __device__ const auto device_obj<RValueMatrix<Derived>>::col(size_t c) const noexcept {
        return ColVector(Base::getConstCastDerived(), 0, getRow(), c);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::rows(size_t fromRow, size_t rowCount) noexcept {
        return BlockType(Base::getDerived(), fromRow, rowCount, 0, getCol());
    }

    template<class Derived>
    __host__ __device__ const auto device_obj<RValueMatrix<Derived>>::rows(size_t fromRow, size_t rowCount) const noexcept {
        return BlockType(Base::getConstCastDerived(), fromRow, rowCount, 0, getCol());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::topRows(size_t to) noexcept {
        return BlockType(Base::getDerived(), 0, to, 0, getCol());
    }

    template<class Derived>
    __host__ __device__ const auto device_obj<RValueMatrix<Derived>>::topRows(size_t to) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, to, 0, getCol());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::bottomRows(size_t from) noexcept {
        return BlockType(Base::getDerived(), from, getRow() - from, 0, getCol());
    }

    template<class Derived>
    __host__ __device__ const auto device_obj<RValueMatrix<Derived>>::bottomRows(size_t from) const noexcept {
        return BlockType(Base::getConstCastDerived(), from, getRow() - from, 0, getCol());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::cols(size_t fromCol, size_t colCount) noexcept {
        return BlockType(Base::getDerived(), 0, getRow(), fromCol, colCount);
    }

    template<class Derived>
    __host__ __device__ const auto device_obj<RValueMatrix<Derived>>::cols(size_t fromCol, size_t colCount) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, getRow(), fromCol, colCount);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::leftCols(size_t to) noexcept {
        return BlockType(Base::getDerived(), 0, getRow(), 0, to);
    }

    template<class Derived>
    __host__ __device__ const auto device_obj<RValueMatrix<Derived>>::leftCols(size_t to) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, getRow(), 0, to);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::rightCols(size_t from) noexcept {
        return BlockType(Base::getDerived(), 0, getRow(), from, getCol() - from);
    }

    template<class Derived>
    __host__ __device__ const auto device_obj<RValueMatrix<Derived>>::rightCols(size_t from) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, getRow(), from, getCol() - from);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::topLeftCorner(size_t toRow, size_t toCol) noexcept {
        return BlockType(Base::getDerived(), 0, toRow, 0, toCol);
    }

    template<class Derived>
    __host__ __device__ const auto device_obj<RValueMatrix<Derived>>::topLeftCorner(size_t toRow, size_t toCol) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, toRow, 0, toCol);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::topLeftCorner(size_t to) noexcept {
        return BlockType(Base::getDerived(), 0, to, 0, to);
    }

    template<class Derived>
    __host__ __device__ const auto device_obj<RValueMatrix<Derived>>::topLeftCorner(size_t to) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, to, 0, to);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::topRightCorner(size_t toRow, size_t fromCol) noexcept {
        return BlockType(Base::getDerived(), 0, toRow, fromCol, getRow() - fromCol);
    }

    template<class Derived>
    __host__ __device__ const auto device_obj<RValueMatrix<Derived>>::topRightCorner(size_t toRow, size_t fromCol) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, toRow, fromCol, getRow() - fromCol);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::bottomLeftCorner(size_t fromRow, size_t toCol) noexcept {
        return BlockType(Base::getDerived(), fromRow, getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    __host__ __device__ const auto device_obj<RValueMatrix<Derived>>::bottomLeftCorner(size_t fromRow, size_t toCol) const noexcept {
        return BlockType(Base::getConstCastDerived(), fromRow, getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::bottomRightCorner(size_t fromRow, size_t fromCol) noexcept {
        return BlockType(Base::getDerived(), fromRow, getRow() - fromRow, fromCol, getCol() - fromCol);
    }

    template<class Derived>
    __host__ __device__ const auto device_obj<RValueMatrix<Derived>>::bottomRightCorner(size_t fromRow, size_t fromCol) const noexcept {
        return BlockType(Base::getConstCastDerived(), fromRow, getRow() - fromRow, fromCol, getCol() - fromCol);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::bottomRightCorner(size_t from) noexcept {
        return BlockType(Base::getDerived(), from, getRow() - from, from, getCol() - from);
    }

    template<class Derived>
    __host__ __device__ const auto device_obj<RValueMatrix<Derived>>::bottomRightCorner(size_t from) const noexcept {
        return BlockType(Base::getConstCastDerived(), from, getRow() - from, from, getCol() - from);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept {
        return BlockType(Base::getDerived(), fromRow, rowCount, fromCol, colCount);
    }

    template<class Derived>
    __host__ __device__ const auto device_obj<RValueMatrix<Derived>>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) const noexcept {
        return BlockType(Base::getConstCastDerived(), fromRow, rowCount, fromCol, colCount);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::diag(this auto&& self) noexcept {
        assert(self.isSquare());
        using Self = decltype(self);
        return device_obj<DiagVectorR<Derived>>(std::forward<Self>(self));
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::triu(this auto&& self) noexcept {
        using Self = decltype(self);
        return device_obj<MatrixTrig<Derived, true, false>>(std::forward<Self>(self));
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::triu_unit(this auto&& self) noexcept {
        using Self = decltype(self);
        return device_obj<MatrixTrig<Derived, true, true>>(std::forward<Self>(self));
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::tril(this auto&& self) noexcept {
        using Self = decltype(self);
        return device_obj<MatrixTrig<Derived, false, false>>(std::forward<Self>(self));
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::tril_unit(this auto&& self) noexcept {
        using Self = decltype(self);
        return device_obj<MatrixTrig<Derived, false, true>>(std::forward<Self>(self));
    }

    template<class Derived>
    __device__ auto device_obj<RValueMatrix<Derived>>::calcFromMajorMinor(size_t major, size_t minor) const -> T {
        return calc(rowFromMajorMinor(major, minor), colFromMajorMinor(major, minor));
    }

    template<class Derived>
    void device_obj<RValueMatrix<Derived>>::reverse(const Matrix auto&, const Matrix auto& grad) const noexcept {
        static_assert(isReverseDiff);
        Base::getDerived().reverse(grad);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::sum_rows() const {
        return device_obj<MatrixSum<Derived, false>>(Base::getDerived());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::sum_cols() const {
        return device_obj<MatrixSum<Derived, true>>(Base::getDerived());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::transpose() const noexcept {
        return device_obj<Transpose<Derived>>(Base::getDerived());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::hermite() const noexcept {
        if constexpr (isComplex)
            return device_obj<Hermite<Derived>>(Base::getDerived());
        else
            return transpose();
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::flatten() const noexcept {
        return device_obj<FlattenR<Derived>>(Base::getDerived());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::reals() const noexcept -> RealsRtnTy {
        return RealsRtnTy(Base::getDerived());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::imags() const noexcept {
        return device_obj<ImagMatrix<Derived>>(Base::getDerived());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::squaredNorms() const noexcept {
        return device_obj<SquaredNormMatrix<Derived>>(Base::getDerived());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::norms() const noexcept {
        return device_obj<NormMatrix<Derived>>(Base::getDerived());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::values() const noexcept -> ValuesRtnTy {
        return ValuesRtnTy(Base::getDerived());
    }

    template<class Derived>
    template<int GradOrder>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::grads() const noexcept {
        return device_obj<GradMatrix<Derived, GradOrder>>(Base::getDerived());
    }

    template<class Derived>
    __host__ __device__ KernelConfig device_obj<RValueMatrix<Derived>>::makeKernelConfig() const noexcept {
        return makeKernelConfig(getMaxMajor(), getMaxMinor());
    }

    template<class Derived>
    __host__ __device__ bool device_obj<RValueMatrix<Derived>>::isSquare() const noexcept {
        return getRow() == getCol();
    }

    template<class Derived>
    __host__ __device__ size_t device_obj<RValueMatrix<Derived>>::rowFromMajorMinor(size_t major, size_t minor) noexcept {
        return MatrixOption::rowFromMajorMinor<device_obj<Derived>>(major, minor);
    }

    template<class Derived>
    __host__ __device__ size_t device_obj<RValueMatrix<Derived>>::colFromMajorMinor(size_t major, size_t minor) noexcept {
        return MatrixOption::colFromMajorMinor<device_obj<Derived>>(major, minor);
    }

    template<class Derived>
    __host__ __device__ KernelConfig device_obj<RValueMatrix<Derived>>::makeKernelConfig(size_t maxMajor, size_t maxMinor) noexcept {
        constexpr size_t MaxThread = MaxThreadsPerBlock;
        const uint32_t numThread = std::min<uint32_t>(maxMinor, MaxThread);
        const uint32_t numBlockX = (maxMinor + numThread - 1) / numThread;
        const uint32_t numBlockY = maxMajor;
        return KernelConfig({numBlockX, numBlockY}, numThread);
    }

    template<class Derived>
    __host__ __device__ void device_obj<RValueMatrix<Derived>>::static_assert_assign(const Matrix auto& source) noexcept {
        static_assert(SizeAtCompile != Dynamic || CUDA<decltype(source)>, "[Error]: Host object cannot be assigned to dynamic device object");
        host_obj::static_assert_assign(source);
    }
}

#include "Trig/MatrixTrig.cuh"
