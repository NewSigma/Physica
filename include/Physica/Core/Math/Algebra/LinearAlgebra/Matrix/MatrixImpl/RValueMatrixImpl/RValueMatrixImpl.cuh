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
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::operator*(this auto&& self, Vector auto&& v) noexcept {
        using Self = decltype(self);
        using V = decltype(v);
        static_assert(is_device_obj<V>::value, "[Error]: host-device mismatch");
        return device_obj<GEMV<remove_device_obj_t<Self>, remove_device_obj_t<V>>>(std::forward<Self>(self), std::forward<V>(v));
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::operator*(this auto&& self, Matrix auto&& m) noexcept {
        using Self = decltype(self);
        using M = decltype(m);
        static_assert(is_device_obj<M>::value, "[Error]: host-device mismatch");
        constexpr bool ColVectorLHS = self.getColAtCompile() == 1;
        constexpr bool RowVectorLHS = self.getRowAtCompile() == 1;
        constexpr bool ColVectorRHS = m.getColAtCompile() == 1;
        constexpr bool RowVectorRHS = m.getRowAtCompile() == 1;
        if constexpr (RowVectorLHS)
            return (std::forward<M>(m).transpose() * std::forward<Self>(self).row(0)).transpose();
        else if constexpr (ColVectorRHS)
            return std::forward<Self>(self) * std::forward<M>(m).col(0);
        else if constexpr (ColVectorLHS || RowVectorRHS) {
            assert(self.getCol() == m.getRow());
            return std::forward<Self>(self).col(0) * std::forward<M>(m);
        }
        else
            return device_obj<GEMM<remove_device_obj_t<Self>, remove_device_obj_t<M>>>(std::forward<Self>(self), std::forward<M>(m));
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::operator-(this auto&& self) noexcept {
        using Self = decltype(self);
        return device_obj<MatrixExpr<ExprID::Minus, remove_device_obj_t<Self>>>(std::forward<Self>(self));
    }

    template<class Derived>
    __host__ __device__ void device_obj<RValueMatrix<Derived>>::assign(Matrix auto&& target) const {
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
            CUDAExecutor::launch<CUDADevAttr::DefaultThreadsPerBlock>(func, target.makeKernelConfig());
        }

        if constexpr (IsDevice())
            Base::getDerived().assign(target, ThreadBlock<1>{});
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
        if constexpr (std::same_as<device_obj<Derived>, std::remove_cvref_t<decltype(source)>>)
            assert(this != &source && "[Error]: Self assign is likely a bug");

        constexpr size_t Row1 = Derived::getRowAtCompile();
        constexpr size_t Row2 = source.getRowAtCompile();
        if constexpr (Row1 == Dynamic || Row2 == Dynamic)
            assert(getRow() == source.getRow() && "[Error]: Dimensions do not match");

        constexpr size_t Col1 = Derived::getColAtCompile();
        constexpr size_t Col2 = source.getColAtCompile();
        if constexpr (Col1 == Dynamic || Col2 == Dynamic)
            assert(getCol() == source.getCol() && "[Error]: Dimensions do not match");
        
        constexpr size_t Size1 = Derived::getSizeAtCompile();
        constexpr size_t Size2 = source.getSizeAtCompile();
        if constexpr (Size1 == Dynamic || Size2 == Dynamic)
            assert(getSize() > 0);
    }

    template<class Derived>
    __device__ void device_obj<RValueMatrix<Derived>>::assign(this const auto& self, Matrix auto&& target, instanceof_x<ThreadBlock> auto block) {
        const size_t maxMinor = target.getMaxMinor();
        for (unsigned int i = block.tid(); i < target.getSize(); i += block.getNumThread()) {
            size_t major = i / maxMinor;
            size_t minor = i % maxMinor;
            const size_t r = target.rowFromMajorMinor(major, minor);
            const size_t c = target.colFromMajorMinor(major, minor);
            target[r, c] = self.calc(r, c);
        }
        block.sync();
    }

    template<class Derived>
    __device__ void device_obj<RValueMatrix<Derived>>::assign_add(this const auto& self, Matrix auto&& target, instanceof_x<ThreadBlock> auto block) {
        const size_t maxMinor = target.getMaxMinor();
        for (unsigned int i = block.tid(); i < target.getSize(); i += block.getNumThread()) {
            size_t major = i / maxMinor;
            size_t minor = i % maxMinor;
            const size_t r = target.rowFromMajorMinor(major, minor);
            const size_t c = target.colFromMajorMinor(major, minor);
            target[r, c] += self.calc(r, c);
        }
        block.sync();
    }

    template<class Derived>
    __host__ __device__ constexpr KernelConfig  device_obj<RValueMatrix<Derived>>::makeKernelConfig() const noexcept {
        constexpr size_t MaxThread = CUDADevAttr::DefaultThreadsPerBlock;
        const uint32_t numThread = std::min<uint32_t>(getMaxMinor(), MaxThread);
        const uint32_t numBlockX = (getMaxMinor() + numThread - 1) / numThread;
        const uint32_t numBlockY = getMaxMajor();
        return KernelConfig ({numBlockX, numBlockY}, numThread);
    }

    template<class Derived>
    __device__ auto device_obj<RValueMatrix<Derived>>::calc(size_t row, size_t col) const {
        return calc(row, col, ThreadBlock<1>{});
    }

    template<class Derived>
    __device__ auto device_obj<RValueMatrix<Derived>>::calc(size_t row, size_t col, instanceof_x<ThreadBlock> auto block) const {
        return Base::getDerived().calc(row, col, block);
    }

    template<class Derived>
    __device__ auto device_obj<RValueMatrix<Derived>>::calc_value(size_t row, size_t col) const {
        return calc_value(row, col, ThreadBlock<1>{});
    }

    template<class Derived>
    __device__ auto device_obj<RValueMatrix<Derived>>::calc_value(size_t row, size_t col, instanceof_x<ThreadBlock> auto block) const {
        return Base::getDerived().values().calc_value(row, col, block);
    }

    template<class Derived>
    __device__ auto device_obj<RValueMatrix<Derived>>::calcFromMajorMinor(this const auto& self, size_t major, size_t minor) -> T {
        return self.calc(rowFromMajorMinor(major, minor), colFromMajorMinor(major, minor));
    }

    template<class Derived>
    void device_obj<RValueMatrix<Derived>>::reverse(const Matrix auto&, const Matrix auto& grad) const noexcept {
        static_assert(isReverseDiff());
        Base::getDerived().reverse(grad);
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
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::row(this auto&& self, size_t r) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<RMatrixBlock<M, 1, Dynamic>>(std::forward<Self>(self), r, 0, self.getCol());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::col(this auto&& self, size_t c) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<RMatrixBlock<M, Dynamic, 1>>(std::forward<Self>(self), 0, self.getRow(), c);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::rows(this auto&& self, size_t fromRow, size_t rowCount) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<RMatrixBlock<M>>(std::forward<Self>(self), fromRow, rowCount, 0, self.getCol());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::topRows(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<RMatrixBlock<M>>(std::forward<Self>(self), 0, to, 0, self.getCol());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::bottomRows(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<RMatrixBlock<M>>(std::forward<Self>(self), from, self.getRow() - from, 0, self.getCol());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::cols(this auto&& self, size_t fromCol, size_t colCount) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<RMatrixBlock<M>>(std::forward<Self>(self), 0, self.getRow(), fromCol, colCount);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::leftCols(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<RMatrixBlock<M>>(std::forward<Self>(self), 0, self.getRow(), 0, to);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::rightCols(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<RMatrixBlock<M>>(std::forward<Self>(self), 0, self.getRow(), from, self.getCol() - from);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::topLeftCorner(this auto&& self, size_t toRow, size_t toCol) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<RMatrixBlock<M>>(std::forward<Self>(self), 0, toRow, 0, toCol);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::topLeftCorner(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<RMatrixBlock<M>>(std::forward<Self>(self), 0, to, 0, to);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::topRightCorner(this auto&& self, size_t toRow, size_t fromCol) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<RMatrixBlock<M>>(std::forward<Self>(self), 0, toRow, fromCol, self.getRow() - fromCol);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::bottomLeftCorner(this auto&& self, size_t fromRow, size_t toCol) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<RMatrixBlock<M>>(std::forward<Self>(self), fromRow, self.getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::bottomRightCorner(this auto&& self, size_t fromRow, size_t fromCol) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<RMatrixBlock<M>>(std::forward<Self>(self), fromRow, self.getRow() - fromRow, fromCol, self.getCol() - fromCol);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::bottomRightCorner(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<RMatrixBlock<M>>(std::forward<Self>(self), from, self.getRow() - from, from, self.getCol() - from);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::block(this auto&& self, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<RMatrixBlock<M>>(std::forward<Self>(self), fromRow, rowCount, fromCol, colCount);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::diag(this auto&& self) noexcept {
        assert(self.isSquare());
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<DiagVectorR<M>>(std::forward<Self>(self));
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::diag(this auto&& self, ssize_t shift) noexcept {
        assert(self.isSquare());
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<MinorDiagR<M>>(std::forward<Self>(self), shift);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::triu(this auto&& self) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<MatrixTrig<M, true, false>>(std::forward<Self>(self));
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::triu_unit(this auto&& self) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<MatrixTrig<M, true, true>>(std::forward<Self>(self));
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::tril(this auto&& self) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<MatrixTrig<M, false, false>>(std::forward<Self>(self));
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::tril_unit(this auto&& self) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<MatrixTrig<M, false, true>>(std::forward<Self>(self));
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
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::inv(this auto&& self) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<Inverse<M>>(std::forward<Self>(self));
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::transpose(this auto&& self) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<Transpose<M>>(std::forward<Self>(self));
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::hermite() const noexcept {
        if constexpr (isComplex())
            return device_obj<Hermite<Derived>>(Base::getDerived());
        else
            return transpose();
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::flatten() const noexcept {
        return device_obj<FlattenR<Derived>>(Base::getDerived());
    }

    template<class Derived>
    __host__ __device__ decltype(auto) device_obj<RValueMatrix<Derived>>::reals(this auto&& self) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        if constexpr (isComplex())
            return device_obj<RealMatrix<M>>(std::forward<Self>(self));
        else
            return std::forward<Self>(self);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::imags(this auto&& self) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<ImagMatrix<M>>(std::forward<Self>(self));
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::squaredNorms(this auto&& self) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<SquaredNormMatrix<M>>(std::forward<Self>(self));
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::norms(this auto&& self) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<NormMatrix<M>>(std::forward<Self>(self));
    }

    template<class Derived>
    __host__ __device__ decltype(auto) device_obj<RValueMatrix<Derived>>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        if constexpr (isDiffable())
            return device_obj<ValueMatrix<M>>(std::forward<Self>(self));
        else
            return std::forward<Self>(self);
    }

    template<class Derived>
    template<int GradOrder>
    __host__ __device__ auto device_obj<RValueMatrix<Derived>>::grads(this auto&& self) noexcept {
        using Self = decltype(self);
        using M = remove_device_obj<Self>::type;
        return device_obj<GradMatrix<M, GradOrder>>(std::forward<Self>(self));
    }

    template<class Derived>
    __host__ __device__ bool device_obj<RValueMatrix<Derived>>::isSquare() const noexcept {
        return getRow() == getCol();
    }

    template<class Derived>
    __host__ __device__ consteval bool device_obj<RValueMatrix<Derived>>::isComplex() noexcept {
        return ScalarType::isComplex();
    }

    template<class Derived>
    __host__ __device__ consteval bool device_obj<RValueMatrix<Derived>>::isDiffable() noexcept {
        return ScalarType::isDiffable();
    }

    template<class Derived>
    __host__ __device__ consteval bool device_obj<RValueMatrix<Derived>>::isForwardDiff() noexcept {
        return ScalarType::isForwardDiff();
    }

    template<class Derived>
    __host__ __device__ consteval bool device_obj<RValueMatrix<Derived>>::isReverseDiff() noexcept {
        return ScalarType::isReverseDiff();
    }

    template<class Derived>
    __host__ __device__ consteval bool device_obj<RValueMatrix<Derived>>::isLValueMatrix() noexcept {
        return false;
    }

    template<class Derived>
    __host__ __device__ consteval bool device_obj<RValueMatrix<Derived>>::isCompact() noexcept {
        return host_obj::isCompact();
    }

    template<class Derived>
    __host__ __device__ consteval bool device_obj<RValueMatrix<Derived>>::isSparse() noexcept {
        return host_obj::isSparse();
    }

    template<class Derived>
    __host__ __device__ consteval bool device_obj<RValueMatrix<Derived>>::isStaticSymm() noexcept {
        return host_obj::isStaticSymm();
    }

    template<class Derived>
    __host__ __device__ consteval bool device_obj<RValueMatrix<Derived>>::isStaticHermite() noexcept {
        return host_obj::isStaticHermite();
    }

    template<class Derived>
    __host__ __device__ consteval bool device_obj<RValueMatrix<Derived>>::isColMatrix() noexcept {
        return host_obj::isColMatrix();
    }

    template<class Derived>
    __host__ __device__ consteval bool device_obj<RValueMatrix<Derived>>::isRowMatrix() noexcept {
        return host_obj::isRowMatrix();
    }

    template<class Derived>
    __host__ __device__ consteval bool device_obj<RValueMatrix<Derived>>::isBothMajor() noexcept {
        return host_obj::isBothMajor();
    }

    template<class Derived>
    __host__ __device__ consteval bool device_obj<RValueMatrix<Derived>>::isFastAssign() noexcept {
        return Derived::isFastAssign();
    }

    template<class Derived>
    __host__ __device__ consteval size_t device_obj<RValueMatrix<Derived>>::getRowAtCompile() noexcept {
        return Derived::getRowAtCompile();
    }

    template<class Derived>
    __host__ __device__ consteval size_t device_obj<RValueMatrix<Derived>>::getColAtCompile() noexcept {
        return Derived::getColAtCompile();
    }

    template<class Derived>
    __host__ __device__ consteval size_t device_obj<RValueMatrix<Derived>>::getSizeAtCompile() noexcept {
        return Derived::getSizeAtCompile();
    }

    template<class Derived>
    __host__ __device__ consteval int device_obj<RValueMatrix<Derived>>::getMajor() noexcept {
        return Derived::getMajor();
    }

    template<class Derived>
    __host__ __device__ size_t device_obj<RValueMatrix<Derived>>::rowFromMajorMinor(size_t major, size_t minor) noexcept {
        return MatrixMajor::rowFromMajorMinor<device_obj<Derived>>(major, minor);
    }

    template<class Derived>
    __host__ __device__ size_t device_obj<RValueMatrix<Derived>>::colFromMajorMinor(size_t major, size_t minor) noexcept {
        return MatrixMajor::colFromMajorMinor<device_obj<Derived>>(major, minor);
    }
    // Redeclare to expose it to base classes
    template<class Derived>
    __host__ __device__ consteval void device_obj<RValueMatrix<Derived>>::static_assert_assign(const Scalar auto& source) noexcept {
        host_obj::static_assert_assign(source);
    }

    template<class Derived>
    __host__ __device__ consteval void device_obj<RValueMatrix<Derived>>::static_assert_assign(const Matrix auto& source) noexcept {
        static_assert(getSizeAtCompile() != Dynamic || DeviceObj<decltype(source)>, "[Error]: Host object cannot be assigned to dynamic device object");
        host_obj::static_assert_assign(source);
    }
}

#include "Trig/MatrixTrig.cuh"
#include "Inverse.cuh"
