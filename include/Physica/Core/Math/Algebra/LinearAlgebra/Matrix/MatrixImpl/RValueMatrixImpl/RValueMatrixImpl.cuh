/*
 * Copyright 2022-2024 Weibo He.
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

#include <Physica/Core/Exception/CUDA/CUDA.cuh>
#include <Physica/Core/Parallel/CUDAContext.cuh>
#include <Physica/PlainStruct.h>
#include "RValueFlatten.cuh"

namespace Physica::Core {
    namespace Internal {
        template<class Derived, Matrix M>
        __global__ void RValueMatrix_assignToKernel(
                Physica::PlainStruct<const device_obj<Derived>> source_, Physica::PlainStruct<device_obj<M>> target_) {
            using ScalarType = typename device_obj<M>::ScalarType;
            const auto& source = source_.getDerived();
            auto& target = target_.getDerived();
            const size_t maxMinor = source.getMaxMinor();
            const size_t major = blockIdx.y;
            const size_t minor = blockIdx.x * blockDim.x + threadIdx.x;
            if (minor < maxMinor) {
                const size_t row = target.rowFromMajorMinor(major, minor);
                const size_t col = target.colFromMajorMinor(major, minor);
                target.refFromMajorMinor(major, minor) = ScalarType(source.calc(row, col));
            }
        }
    }

    template<class Derived>
    template<Matrix M>
    __host__ __device__ void device_obj<RValueMatrix<Derived>>::assignTo(device_obj<LValueMatrix<M>>& target) const {
        [[maybe_unused]] const auto kernel = Internal::RValueMatrix_assignToKernel<Derived, M>;
        const size_t maxMajor = target.getMaxMajor();
        const size_t maxMinor = target.getMaxMinor();
        if constexpr (IsHost()) {
            const unsigned int numThread = std::min(maxMinor, MaxThreadPerBlock);
            const unsigned int numBlockX = (maxMinor + numThread - 1) / numThread;
            const unsigned int numBlockY = maxMajor;
            kernel<<<{numBlockX, numBlockY}, numThread, 0, CUDAContext::getInstance()>>>(asStruct(Base::getDerived()), asStruct(target.getDerived()));
            check(cudaGetLastError());
        }
        else {
            using OtherScalar = typename M::ScalarType;
            for (size_t major = 0; major < maxMajor; ++major) {
                for (size_t minor = 0; minor < maxMinor; ++minor) {
                    const size_t r = target.rowFromMajorMinor(major, minor);
                    const size_t c = target.colFromMajorMinor(major, minor);
                    target.refFromMajorMinor(major, minor) = OtherScalar(calc(r, c));
                }
            }
        }
    }

    template<class Derived>
    __host__ __device__ inline typename device_obj<RValueMatrix<Derived>>::RowVector
    device_obj<RValueMatrix<Derived>>::row(size_t r) {
        return {Base::getDerived(), r, 0, getCol()};
    }

    template<class Derived>
    __host__ __device__ inline const typename device_obj<RValueMatrix<Derived>>::RowVector
    device_obj<RValueMatrix<Derived>>::row(size_t r) const {
        return {Base::getConstCastDerived(), r, 0, getCol()};
    }

    template<class Derived>
    __host__ __device__ inline typename device_obj<RValueMatrix<Derived>>::ColVector
    device_obj<RValueMatrix<Derived>>::col(size_t c) {
        return {Base::getDerived(), 0, getRow(), c};
    }

    template<class Derived>
    __host__ __device__ inline const typename device_obj<RValueMatrix<Derived>>::ColVector
    device_obj<RValueMatrix<Derived>>::col(size_t c) const {
        return {Base::getConstCastDerived(), 0, getRow(), c};
    }

    template<class Derived>
    __host__ __device__ inline device_obj<RValueMatrix<Derived>>::BlockType
    device_obj<RValueMatrix<Derived>>::rows(size_t fromRow, size_t rowCount) {
        return {Base::getDerived(), fromRow, rowCount, 0, getCol()};
    }

    template<class Derived>
    __host__ __device__ inline const device_obj<RValueMatrix<Derived>>::BlockType
    device_obj<RValueMatrix<Derived>>::rows(size_t fromRow, size_t rowCount) const {
        return {Base::getConstCastDerived(), fromRow, rowCount, 0, getCol()};
    }

    template<class Derived>
    __host__ __device__ inline device_obj<RValueMatrix<Derived>>::BlockType
    device_obj<RValueMatrix<Derived>>::topRows(size_t to) {
        return {Base::getDerived(), 0, to, 0, getCol()};
    }

    template<class Derived>
    __host__ __device__ inline const device_obj<RValueMatrix<Derived>>::BlockType
    device_obj<RValueMatrix<Derived>>::topRows(size_t to) const {
        return {Base::getConstCastDerived(), 0, to, 0, getCol()};
    }

    template<class Derived>
    __host__ __device__ inline device_obj<RValueMatrix<Derived>>::BlockType
    device_obj<RValueMatrix<Derived>>::bottomRows(size_t from) {
        return {Base::getDerived(), from, getRow() - from, 0, getCol()};
    }

    template<class Derived>
    __host__ __device__ inline const device_obj<RValueMatrix<Derived>>::BlockType
    device_obj<RValueMatrix<Derived>>::bottomRows(size_t from) const {
        return {Base::getConstCastDerived(), from, getRow() - from, 0, getCol()};
    }

    template<class Derived>
    __host__ __device__ inline device_obj<RValueMatrix<Derived>>::BlockType
    device_obj<RValueMatrix<Derived>>::cols(size_t fromCol, size_t colCount) {
        return {Base::getDerived(), 0, getRow(), fromCol, colCount};
    }

    template<class Derived>
    __host__ __device__ inline const device_obj<RValueMatrix<Derived>>::BlockType
    device_obj<RValueMatrix<Derived>>::cols(size_t fromCol, size_t colCount) const {
        return {Base::getConstCastDerived(), 0, getRow(), fromCol, colCount};
    }

    template<class Derived>
    __host__ __device__ inline device_obj<RValueMatrix<Derived>>::BlockType
    device_obj<RValueMatrix<Derived>>::leftCols(size_t to) {
        return {Base::getDerived(), 0, getRow(), 0, to};
    }

    template<class Derived>
    __host__ __device__ inline const device_obj<RValueMatrix<Derived>>::BlockType
    device_obj<RValueMatrix<Derived>>::leftCols(size_t to) const {
        return {Base::getConstCastDerived(), 0, getRow(), 0, to};
    }

    template<class Derived>
    __host__ __device__ inline device_obj<RValueMatrix<Derived>>::BlockType
    device_obj<RValueMatrix<Derived>>::rightCols(size_t from) {
        return {Base::getDerived(), 0, getRow(), from, getCol() - from};
    }

    template<class Derived>
    __host__ __device__ inline const device_obj<RValueMatrix<Derived>>::BlockType
    device_obj<RValueMatrix<Derived>>::rightCols(size_t from) const {
        return {Base::getConstCastDerived(), 0, getRow(), from, getCol() - from};
    }

    template<class Derived>
    __host__ __device__ inline device_obj<RValueMatrix<Derived>>::BlockType
    device_obj<RValueMatrix<Derived>>::topLeftCorner(size_t toRow, size_t toCol) {
        return {Base::getDerived(), 0, toRow, 0, toCol};
    }

    template<class Derived>
    __host__ __device__ inline const device_obj<RValueMatrix<Derived>>::BlockType
    device_obj<RValueMatrix<Derived>>::topLeftCorner(size_t toRow, size_t toCol) const {
        return {Base::getConstCastDerived(), 0, toRow, 0, toCol};
    }

    template<class Derived>
    __host__ __device__ inline device_obj<RValueMatrix<Derived>>::BlockType
    device_obj<RValueMatrix<Derived>>::topLeftCorner(size_t to) {
        return {Base::getDerived(), 0, to, 0, to};
    }

    template<class Derived>
    __host__ __device__ inline const device_obj<RValueMatrix<Derived>>::BlockType
    device_obj<RValueMatrix<Derived>>::topLeftCorner(size_t to) const {
        return {Base::getConstCastDerived(), 0, to, 0, to};
    }

    template<class Derived>
    __host__ __device__ inline device_obj<RValueMatrix<Derived>>::BlockType
    device_obj<RValueMatrix<Derived>>::topRightCorner(size_t toRow, size_t fromCol) {
        return {Base::getDerived(), 0, toRow, fromCol, getRow() - fromCol};
    }

    template<class Derived>
    __host__ __device__ inline const device_obj<RValueMatrix<Derived>>::BlockType
    device_obj<RValueMatrix<Derived>>::topRightCorner(size_t toRow, size_t fromCol) const {
        return {Base::getConstCastDerived(), 0, toRow, fromCol, getRow() - fromCol};
    }

    template<class Derived>
    __host__ __device__ inline device_obj<RValueMatrix<Derived>>::BlockType
    device_obj<RValueMatrix<Derived>>::bottomLeftCorner(size_t fromRow, size_t toCol) {
        return {Base::getDerived(), fromRow, getRow() - fromRow, 0, toCol};
    }

    template<class Derived>
    __host__ __device__ inline const device_obj<RValueMatrix<Derived>>::BlockType
    device_obj<RValueMatrix<Derived>>::bottomLeftCorner(size_t fromRow, size_t toCol) const {
        return {Base::getConstCastDerived(), fromRow, getRow() - fromRow, 0, toCol};
    }

    template<class Derived>
    __host__ __device__ inline device_obj<RValueMatrix<Derived>>::BlockType
    device_obj<RValueMatrix<Derived>>::bottomRightCorner(size_t fromRow, size_t fromCol) {
        return {Base::getDerived(), fromRow, getRow() - fromRow, fromCol, getCol() - fromCol};
    }

    template<class Derived>
    __host__ __device__ inline const device_obj<RValueMatrix<Derived>>::BlockType
    device_obj<RValueMatrix<Derived>>::bottomRightCorner(size_t fromRow, size_t fromCol) const {
        return {Base::getConstCastDerived(), fromRow, getRow() - fromRow, fromCol, getCol() - fromCol};
    }

    template<class Derived>
    __host__ __device__ inline device_obj<RValueMatrix<Derived>>::BlockType
    device_obj<RValueMatrix<Derived>>::bottomRightCorner(size_t from) {
        return {Base::getDerived(), from, getRow() - from, from, getCol() - from};
    }

    template<class Derived>
    __host__ __device__ inline const device_obj<RValueMatrix<Derived>>::BlockType
    device_obj<RValueMatrix<Derived>>::bottomRightCorner(size_t from) const {
        return {Base::getConstCastDerived(), from, getRow() - from, from, getCol() - from};
    }

    template<class Derived>
    __host__ __device__ inline device_obj<RValueMatrix<Derived>>::BlockType
    device_obj<RValueMatrix<Derived>>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) {
        return {Base::getDerived(), fromRow, rowCount, fromCol, colCount};
    }

    template<class Derived>
    __host__ __device__ inline const device_obj<RValueMatrix<Derived>>::BlockType
    device_obj<RValueMatrix<Derived>>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) const {
        return {Base::getConstCastDerived(), fromRow, rowCount, fromCol, colCount};
    }

    template<class Derived>
    inline typename device_obj<RValueMatrix<Derived>>::ScalarType
    __device__ device_obj<RValueMatrix<Derived>>::calcFromMajorMinor(size_t major, size_t minor) const {
        return calc(rowFromMajorMinor(major, minor), colFromMajorMinor(major, minor));
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::transpose() const noexcept {
        return device_obj<Transpose<Derived>>(Base::getDerived());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::flatten() const noexcept {
        return device_obj<RValueFlatten<Derived>>(Base::getDerived());
    }
}
