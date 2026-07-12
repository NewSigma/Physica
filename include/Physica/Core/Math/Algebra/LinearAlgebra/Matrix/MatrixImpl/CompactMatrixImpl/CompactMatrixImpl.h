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

#include "../CompactMatrix.h"

namespace Physica {
    template<class Derived>
    template<ExecutePolicy P>
    void CompactMatrix<Derived>::assign(Matrix auto&& m) const noexcept {
        assign_base<P>(m);
    }

    template<class Derived>
    template<ExecutePolicy P>
    void CompactMatrix<Derived>::assign_base(Matrix auto&& __restrict target) const __restrict noexcept {
        using M = decltype(target);
        const auto& self = Base::getDerived();
        target.assert_assign(self);
        if constexpr (MatrixMajor::isSameMajor<Derived, decltype(target)>()) {
            if constexpr (std::remove_cvref_t<M>::isCompact())
                self.flatten().template assign<P>(target.flatten());
            else
                Base::template assign<P>(target);
        }
        else {
            constexpr size_t BlockingL1 = Base::calcBlockingSize(HostDevAttr::CacheSizeL1D);
            size_t r0 = self.getRow();
            size_t c0 = self.getCol();
            for (size_t r = 0; r < r0; r += BlockingL1) {
                size_t numR = std::min(BlockingL1, r0 - r);
                for (size_t c = 0; c < c0; c += BlockingL1) {
                    size_t numC = std::min(BlockingL1, c0 - c);
                    self.block(r, numR, c, numC).assign(target.block(r, numR, c, numC));
                }
            }
        }
    }

    template<class Derived>
    auto CompactMatrix<Derived>::row(this auto&& self, size_t r) noexcept {
        using Self = decltype(self);
        constexpr size_t ColAtCompile = self.getColAtCompile();
        constexpr bool IsMat1x1 = ColAtCompile == 1;
        if constexpr (IsMat1x1)
            return CompactMatrixBlock<Self, 1, 1>(std::forward<Self>(self), r, 0);
        else {
            if constexpr (Derived::isRowMatrix())
                return CompactMatrixBlock<Self, 1, ColAtCompile>(std::forward<Self>(self), r, 0, self.getCol());
            else
                return LMatrixBlock<Self, 1, ColAtCompile>(std::forward<Self>(self), r, 0, self.getCol());
        }
    }

    template<class Derived>
    auto CompactMatrix<Derived>::col(this auto&& self, size_t c) noexcept {
        using Self = decltype(self);
        constexpr size_t RowAtCompile = self.getRowAtCompile();
        constexpr bool IsMat1x1 = RowAtCompile == 1;
        if constexpr (IsMat1x1)
            return CompactMatrixBlock<Self, 1, 1>(std::forward<Self>(self), 0, c);
        else {
            if constexpr (Derived::isColMatrix())
                return CompactMatrixBlock<Self, RowAtCompile, 1>(std::forward<Self>(self), 0, self.getRow(), c);
            else
                return LMatrixBlock<Self, RowAtCompile, 1>(std::forward<Self>(self), 0, self.getRow(), c);
        }
    }

    template<class Derived>
    template<size_t Row>
    auto CompactMatrix<Derived>::rows(this auto&& self, size_t fromRow, size_t rowCount) noexcept {
        using Self = decltype(self);
        return CompactMatrixBlock<Self, Row, Derived::getColAtCompile()>(std::forward<Self>(self), fromRow, rowCount, 0, self.getCol());
    }

    template<class Derived>
    template<size_t Row>
    auto CompactMatrix<Derived>::topRows(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        return CompactMatrixBlock<Self, Row, Derived::getColAtCompile()>(std::forward<Self>(self), 0, to, 0, self.getCol());
    }

    template<class Derived>
    template<size_t Row>
    auto CompactMatrix<Derived>::bottomRows(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        return CompactMatrixBlock<Self, Row, Derived::getColAtCompile()>(std::forward<Self>(self), from, self.getRow() - from, 0, self.getCol());
    }

    template<class Derived>
    template<size_t Col>
    auto CompactMatrix<Derived>::cols(this auto&& self, size_t fromCol, size_t colCount) noexcept {
        using Self = decltype(self);
        return CompactMatrixBlock<Self, self.getRowAtCompile(), Col>(std::forward<Self>(self), 0, self.getRow(), fromCol, colCount);
    }

    template<class Derived>
    template<size_t Col>
    auto CompactMatrix<Derived>::leftCols(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        return CompactMatrixBlock<Self, self.getRowAtCompile(), Col>(std::forward<Self>(self), 0, self.getRow(), 0, to);
    }

    template<class Derived>
    template<size_t Col>
    auto CompactMatrix<Derived>::rightCols(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        return CompactMatrixBlock<Self, self.getRowAtCompile(), Col>(std::forward<Self>(self), 0, self.getRow(), from, self.getCol() - from);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto CompactMatrix<Derived>::topLeftCorner(this auto&& self, size_t toRow, size_t toCol) noexcept {
        using Self = decltype(self);
        return CompactMatrixBlock<Self, Row, Col>(std::forward<Self>(self), 0, toRow, 0, toCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto CompactMatrix<Derived>::topLeftCorner(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        return CompactMatrixBlock<Self, Row, Col>(std::forward<Self>(self), 0, to, 0, to);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto CompactMatrix<Derived>::topRightCorner(this auto&& self, size_t toRow, size_t fromCol) noexcept {
        using Self = decltype(self);
        return CompactMatrixBlock<Self, Row, Col>(std::forward<Self>(self), 0, toRow, fromCol, self.getRow() - fromCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto CompactMatrix<Derived>::bottomLeftCorner(this auto&& self, size_t fromRow, size_t toCol) noexcept {
        using Self = decltype(self);
        return CompactMatrixBlock<Self, Row, Col>(std::forward<Self>(self), fromRow, self.getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto CompactMatrix<Derived>::bottomRightCorner(this auto&& self, size_t fromRow, size_t fromCol) noexcept {
        using Self = decltype(self);
        return CompactMatrixBlock<Self, Row, Col>(std::forward<Self>(self), fromRow, self.getRow() - fromRow, fromCol, self.getCol() - fromCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto CompactMatrix<Derived>::bottomRightCorner(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        return CompactMatrixBlock<Self, Row, Col>(std::forward<Self>(self), from, self.getRow() - from, from, self.getCol() - from);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto CompactMatrix<Derived>::block(this auto&& self, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept {
        using Self = decltype(self);
        return CompactMatrixBlock<Self, Row, Col>(std::forward<Self>(self), fromRow, rowCount, fromCol, colCount);
    }

    template<class Derived>
    void CompactMatrix<Derived>::zeros() noexcept {
        Base::getDerived().flatten().zeros();
    }

    template<class Derived>
    void CompactMatrix<Derived>::read(const auto& obj) noexcept {
        using O = decltype(obj);
        if constexpr (Vector<O>)
            Base::getDerived().flatten().read(obj);
        else {
            static_assert(Matrix<O>, "[Error]: Unexpected type");
            read(obj.flatten());
        }
    }

    template<class Derived>
    template<RNG R>
    void CompactMatrix<Derived>::random_uniform() {
        Base::getDerived().flatten().template random_uniform<R>();
    }

    template<class Derived>
    template<RNG R>
    void CompactMatrix<Derived>::random_normal() {
        Base::getDerived().flatten().template random_normal<R>();
    }

#ifdef PHYSICA_HDF5
    template<class Derived>
    const H5DataSet<2> CompactMatrix<Derived>::read(const H5Loc& loc, const char* name) {
        const auto dataset = loc.openDataSet<2>(name);
        const size_t maxMajor = dataset.getSize(0);
        const size_t maxMinor = dataset.getSize(1);
        Base::resize(Base::rowFromMajorMinor(maxMajor, maxMinor), Base::colFromMajorMinor(maxMajor, maxMinor));

        const auto memSpace = H5DataSpace<1>({maxMinor});
        auto fileSpace = H5DataSpace<2>({Base::getMaxMajor(), Base::getMaxMinor()});
        for (size_t major = 0; major < maxMajor; ++major) {
            fileSpace.selectHyperslab(H5S_SELECT_SET, {1, maxMinor}, {major, 0});
            if constexpr (Derived::isColMatrix())
                dataset.read(Base::getDerived().col(major).data(), ScalarType::dtype_hdf5(), memSpace, fileSpace);
            else
                dataset.read(Base::getDerived().row(major).data(), ScalarType::dtype_hdf5(), memSpace, fileSpace);
        }
        return dataset;
    }

    template<class Derived>
    H5DataSet<2> CompactMatrix<Derived>::write(H5Loc& loc, const char* name) const {
        const size_t maxMajor = Base::getMaxMajor();
        const size_t maxMinor = Base::getMaxMinor();
        const auto space = H5DataSpace<2>({maxMajor, maxMinor});
        H5DataSet<2> dataset;
        if (loc.exists(name))
            dataset = loc.openDataSet<2>(name);
        else
            dataset = loc.createDataSet<2>(name, ScalarType::dtype_hdf5(), space);

        const auto memSpace = H5DataSpace<1>({maxMinor});
        auto fileSpace = H5DataSpace<2>({Base::getMaxMajor(), Base::getMaxMinor()});
        for (size_t major = 0; major < maxMajor; ++major) {
            fileSpace.selectHyperslab(H5S_SELECT_SET, {1, maxMinor}, {major, 0});
            if constexpr (Derived::isColMatrix())
                dataset.write(Base::getDerived().col(major).data(), ScalarType::dtype_hdf5(), memSpace, fileSpace);
            else
                dataset.write(Base::getDerived().row(major).data(), ScalarType::dtype_hdf5(), memSpace, fileSpace);
        }
        return std::cref(dataset);
    }
#endif

    template<class Derived>
    auto CompactMatrix<Derived>::data() noexcept {
        return Base::getDerived().data();
    }

    template<class Derived>
    auto CompactMatrix<Derived>::data() const noexcept {
        return Base::getDerived().data();
    }

    template<class Derived>
    auto CompactMatrix<Derived>::data_handle(this auto&& self) noexcept {
        return self.data();
    }

    template<class Derived>
    constexpr size_t CompactMatrix<Derived>::getRowStride() const noexcept {
        return Derived::isColMatrix() ? 1 : Base::getCol();
    }

    template<class Derived>
    constexpr size_t CompactMatrix<Derived>::getColStride() const noexcept {
        return Derived::isRowMatrix() ? 1 : Base::getRow();
    }

    template<class Derived>
    __host__ __device__ consteval size_t CompactMatrix<Derived>::getRowStrideAtCompile() noexcept {
        return Derived::isColMatrix() ? 1 : Derived::getColAtCompile();
    }

    template<class Derived>
    __host__ __device__ consteval size_t CompactMatrix<Derived>::getColStrideAtCompile() noexcept {
        return Derived::isRowMatrix() ? 1 : Derived::getRowAtCompile();
    }
}
