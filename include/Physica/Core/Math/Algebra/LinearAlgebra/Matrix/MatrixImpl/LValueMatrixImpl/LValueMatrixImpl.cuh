/*
 * Copyright 2024-2025 Weibo He.
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
    template<Matrix M>
    __host__ __device__ device_obj<Derived>& device_obj<LValueMatrix<Derived>>::operator=(const M& m) requires(CUDA<M>) {
        static_assert(RowAtCompile == Dynamic || M::RowAtCompile == Dynamic || RowAtCompile == M::RowAtCompile, "[Error]: Incompatible row number");
        static_assert(ColAtCompile == Dynamic || M::ColAtCompile == Dynamic || ColAtCompile == M::ColAtCompile, "[Error]: Incompatible col number");
        auto& target = Base::getDerived();
        target.resize(m.getRow(), m.getCol());
        m.assign(target);
        return target;
    }

    template<class Derived>
    template<Scalar T>
    __host__ __device__ device_obj<Derived>& device_obj<LValueMatrix<Derived>>::operator=(const T& x) {
        if (IsHost()) {
            CUDAExecutor::launch([m_ = asStruct(Base::getDerived()), x] __device__() mutable {
                auto& m = m_.getDerived();
                const size_t maxMinor = m.getMaxMinor();
                const size_t major = blockIdx.y;
                const size_t minor = blockIdx.x * blockDim.x + threadIdx.x;
                if (minor < maxMinor)
                    m.refFromMajorMinor(major, minor) = x;
            }, Base::makeKernelConfig());
        }
        else {
            for (size_t i = 0; i < Base::getMaxMajor(); ++i)
                for (size_t j = 0; j < Base::getMaxMinor(); ++j)
                    refFromMajorMinor(i, j) = x;
        }
        return Base::getDerived();
    }

    template<class Derived>
    __device__ auto device_obj<LValueMatrix<Derived>>::operator()(size_t row, size_t col) -> RefTy {
        return *data_ptr(row, col);
    }

    template<class Derived>
    __device__ auto device_obj<LValueMatrix<Derived>>::operator()(size_t row, size_t col) const -> ConstRefTy {
        return *data_ptr(row, col);
    }

    template<class Derived>
    template<Matrix M>
    void device_obj<LValueMatrix<Derived>>::reverse(const M& grad) const noexcept requires(isReverseDiff) {
        static_assert(std::same_as<typename ScalarType::GradType, typename M::ScalarType>, "[Error]: Inconsistent ScalarType");
        assert(Base::getRow() == grad.getRow());
        assert(Base::getCol() == grad.getCol());
        grad.assign_add(Base::getConstCastDerived().grads());
    }

    template<class Derived>
    __host__ __device__ inline auto device_obj<LValueMatrix<Derived>>::row(size_t r) noexcept {
        return RowVector(Base::getDerived(), r, 0, getCol());
    }

    template<class Derived>
    __host__ __device__ inline const auto device_obj<LValueMatrix<Derived>>::row(size_t r) const noexcept {
        return RowVector(Base::getConstCastDerived(), r, 0, getCol());
    }

    template<class Derived>
    __host__ __device__ inline auto device_obj<LValueMatrix<Derived>>::col(size_t c) noexcept {
        return ColVector(Base::getDerived(), 0, getRow(), c);
    }

    template<class Derived>
    __host__ __device__ inline const auto device_obj<LValueMatrix<Derived>>::col(size_t c) const noexcept {
        return ColVector(Base::getConstCastDerived(), 0, getRow(), c);
    }

    template<class Derived>
    __host__ __device__ inline auto device_obj<LValueMatrix<Derived>>::rows(size_t fromRow, size_t rowCount) noexcept {
        return BlockType(Base::getDerived(), fromRow, rowCount, 0, getCol());
    }

    template<class Derived>
    __host__ __device__ inline const auto device_obj<LValueMatrix<Derived>>::rows(size_t fromRow, size_t rowCount) const noexcept {
        return BlockType(Base::getConstCastDerived(), fromRow, rowCount, 0, getCol());
    }

    template<class Derived>
    __host__ __device__ inline auto device_obj<LValueMatrix<Derived>>::topRows(size_t to) noexcept {
        return BlockType(Base::getDerived(), 0, to, 0, getCol());
    }

    template<class Derived>
    __host__ __device__ inline const auto device_obj<LValueMatrix<Derived>>::topRows(size_t to) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, to, 0, getCol());
    }

    template<class Derived>
    __host__ __device__ inline auto device_obj<LValueMatrix<Derived>>::bottomRows(size_t from) noexcept {
        return BlockType(Base::getDerived(), from, getRow() - from, 0, getCol());
    }

    template<class Derived>
    __host__ __device__ inline const auto device_obj<LValueMatrix<Derived>>::bottomRows(size_t from) const noexcept {
        return BlockType(Base::getConstCastDerived(), from, getRow() - from, 0, getCol());
    }

    template<class Derived>
    __host__ __device__ inline auto device_obj<LValueMatrix<Derived>>::cols(size_t fromCol, size_t colCount) noexcept {
        return BlockType(Base::getDerived(), 0, getRow(), fromCol, colCount);
    }

    template<class Derived>
    __host__ __device__ inline const auto device_obj<LValueMatrix<Derived>>::cols(size_t fromCol, size_t colCount) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, getRow(), fromCol, colCount);
    }

    template<class Derived>
    __host__ __device__ inline auto device_obj<LValueMatrix<Derived>>::leftCols(size_t to) noexcept {
        return BlockType(Base::getDerived(), 0, getRow(), 0, to);
    }

    template<class Derived>
    __host__ __device__ inline const auto device_obj<LValueMatrix<Derived>>::leftCols(size_t to) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, getRow(), 0, to);
    }

    template<class Derived>
    __host__ __device__ inline auto device_obj<LValueMatrix<Derived>>::rightCols(size_t from) noexcept {
        return BlockType(Base::getDerived(), 0, getRow(), from, getCol() - from);
    }

    template<class Derived>
    __host__ __device__ inline const auto device_obj<LValueMatrix<Derived>>::rightCols(size_t from) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, getRow(), from, getCol() - from);
    }

    template<class Derived>
    __host__ __device__ inline auto device_obj<LValueMatrix<Derived>>::topLeftCorner(size_t toRow, size_t toCol) noexcept {
        return BlockType(Base::getDerived(), 0, toRow, 0, toCol);
    }

    template<class Derived>
    __host__ __device__ inline const auto device_obj<LValueMatrix<Derived>>::topLeftCorner(size_t toRow, size_t toCol) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, toRow, 0, toCol);
    }

    template<class Derived>
    __host__ __device__ inline auto device_obj<LValueMatrix<Derived>>::topLeftCorner(size_t to) noexcept {
        return BlockType(Base::getDerived(), 0, to, 0, to);
    }

    template<class Derived>
    __host__ __device__ inline const auto device_obj<LValueMatrix<Derived>>::topLeftCorner(size_t to) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, to, 0, to);
    }

    template<class Derived>
    __host__ __device__ inline auto device_obj<LValueMatrix<Derived>>::topRightCorner(size_t toRow, size_t fromCol) noexcept {
        return BlockType(Base::getDerived(), 0, toRow, fromCol, getRow() - fromCol);
    }

    template<class Derived>
    __host__ __device__ inline const auto device_obj<LValueMatrix<Derived>>::topRightCorner(size_t toRow, size_t fromCol) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, toRow, fromCol, getRow() - fromCol);
    }

    template<class Derived>
    __host__ __device__ inline auto device_obj<LValueMatrix<Derived>>::bottomLeftCorner(size_t fromRow, size_t toCol) noexcept {
        return BlockType(Base::getDerived(), fromRow, getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    __host__ __device__ inline const auto device_obj<LValueMatrix<Derived>>::bottomLeftCorner(size_t fromRow, size_t toCol) const noexcept {
        return BlockType(Base::getConstCastDerived(), fromRow, getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    __host__ __device__ inline auto device_obj<LValueMatrix<Derived>>::bottomRightCorner(size_t fromRow, size_t fromCol) noexcept {
        return BlockType(Base::getDerived(), fromRow, getRow() - fromRow, fromCol, getCol() - fromCol);
    }

    template<class Derived>
    __host__ __device__ inline const auto device_obj<LValueMatrix<Derived>>::bottomRightCorner(size_t fromRow, size_t fromCol) const noexcept {
        return BlockType(Base::getConstCastDerived(), fromRow, getRow() - fromRow, fromCol, getCol() - fromCol);
    }

    template<class Derived>
    __host__ __device__ inline auto device_obj<LValueMatrix<Derived>>::bottomRightCorner(size_t from) noexcept {
        return BlockType(Base::getDerived(), from, getRow() - from, from, getCol() - from);
    }

    template<class Derived>
    __host__ __device__ inline const auto device_obj<LValueMatrix<Derived>>::bottomRightCorner(size_t from) const noexcept {
        return BlockType(Base::getConstCastDerived(), from, getRow() - from, from, getCol() - from);
    }

    template<class Derived>
    __host__ __device__ inline auto device_obj<LValueMatrix<Derived>>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept {
        return BlockType(Base::getDerived(), fromRow, rowCount, fromCol, colCount);
    }

    template<class Derived>
    __host__ __device__ inline const auto device_obj<LValueMatrix<Derived>>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) const noexcept {
        return BlockType(Base::getConstCastDerived(), fromRow, rowCount, fromCol, colCount);
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
    __host__ __device__ auto device_obj<LValueMatrix<Derived>>::data_ptr(size_t row, size_t col) -> PtrTy {
        assert(row < Base::getRow());
        assert(col < Base::getCol());
        return Base::getDerived().data_ptr(row, col);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<LValueMatrix<Derived>>::data_ptr(size_t row, size_t col) const -> ConstPtrTy {
        return const_cast<This&>(*this).data_ptr(row, col);
    }

    template<class Derived>
    __device__ inline auto device_obj<LValueMatrix<Derived>>::refFromMajorMinor(size_t major, size_t minor) -> RefTy {
        assert(major < Base::getDerived().getMaxMajor());
        assert(minor < Base::getDerived().getMaxMinor());
        const size_t r = MatrixOption::rowFromMajorMinor<Derived>(major, minor);
        const size_t c = MatrixOption::colFromMajorMinor<Derived>(major, minor);
        return Base::getDerived()(r, c);
    }

    template<class Derived>
    __device__ inline auto device_obj<LValueMatrix<Derived>>::refFromMajorMinor(size_t major, size_t minor) const -> ConstRefTy {
        return const_cast<This&>(*this).refFromMajorMinor(major, minor);
    }
}
