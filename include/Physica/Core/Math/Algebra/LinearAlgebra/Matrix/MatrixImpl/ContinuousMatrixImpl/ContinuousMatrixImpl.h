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

#include "../ContinuousMatrix.h"

namespace Physica {
    template<class Derived>
    template<ExecutePolicy P>
    void ContinuousMatrix<Derived>::assign(Matrix auto&& m) const noexcept {
        using M = decltype(m);
        if constexpr (is_continuous<M>::value && MatrixOption::isSameMajor<Derived, M>())
            Base::getDerived().flatten().template assign<P>(m.flatten());
        else
            Base::template assign<P>(m);
    }

    template<class Derived>
    auto ContinuousMatrix<Derived>::row(this auto&& self, size_t r) noexcept {
        using Self = decltype(self);
        const bool IsMat1x1 = Base::ColAtCompile == 1;
        if constexpr (IsMat1x1)
            return ContinuousMatrixBlock<Self, 1, 1>(std::forward<Self>(self), r, 0);
        else {
            if constexpr (isRowMatrix)
                return ContinuousMatrixBlock<Self, 1, ColAtCompile>(std::forward<Self>(self), r, 0, self.getCol());
            else
                return LMatrixBlock<Self, 1, ColAtCompile>(std::forward<Self>(self), r, 0, self.getCol());
        }
    }

    template<class Derived>
    auto ContinuousMatrix<Derived>::col(this auto&& self, size_t c) noexcept {
        using Self = decltype(self);
        const bool IsMat1x1 = Base::RowAtCompile == 1;
        if constexpr (IsMat1x1)
            return ContinuousMatrixBlock<Self, 1, 1>(std::forward<Self>(self), 0, c);
        else {
            if constexpr (isColMatrix)
                return ContinuousMatrixBlock<Self, RowAtCompile, 1>(std::forward<Self>(self), 0, self.getRow(), c);
            else
                return LMatrixBlock<Self, RowAtCompile, 1>(std::forward<Self>(self), 0, self.getRow(), c);
        }
    }

    template<class Derived>
    template<size_t Row>
    auto ContinuousMatrix<Derived>::rows(this auto&& self, size_t fromRow, size_t rowCount) noexcept {
        using Self = decltype(self);
        return ContinuousMatrixBlock<Self, Row, ColAtCompile>(std::forward<Self>(self), fromRow, rowCount, 0, self.getCol());
    }

    template<class Derived>
    template<size_t Row>
    auto ContinuousMatrix<Derived>::topRows(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        return ContinuousMatrixBlock<Self, Row, ColAtCompile>(std::forward<Self>(self), 0, to, 0, self.getCol());
    }

    template<class Derived>
    template<size_t Row>
    auto ContinuousMatrix<Derived>::bottomRows(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        return ContinuousMatrixBlock<Self, Row, ColAtCompile>(std::forward<Self>(self), from, self.getRow() - from, 0, self.getCol());
    }

    template<class Derived>
    template<size_t Col>
    auto ContinuousMatrix<Derived>::cols(this auto&& self, size_t fromCol, size_t colCount) noexcept {
        using Self = decltype(self);
        return ContinuousMatrixBlock<Self, RowAtCompile, Col>(std::forward<Self>(self), 0, self.getRow(), fromCol, colCount);
    }

    template<class Derived>
    template<size_t Col>
    auto ContinuousMatrix<Derived>::leftCols(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        return ContinuousMatrixBlock<Self, RowAtCompile, Col>(std::forward<Self>(self), 0, self.getRow(), 0, to);
    }

    template<class Derived>
    template<size_t Col>
    auto ContinuousMatrix<Derived>::rightCols(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        return ContinuousMatrixBlock<Self, RowAtCompile, Col>(std::forward<Self>(self), 0, self.getRow(), from, self.getCol() - from);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto ContinuousMatrix<Derived>::topLeftCorner(this auto&& self, size_t toRow, size_t toCol) noexcept {
        using Self = decltype(self);
        return ContinuousMatrixBlock<Self, Row, Col>(std::forward<Self>(self), 0, toRow, 0, toCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto ContinuousMatrix<Derived>::topLeftCorner(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        return ContinuousMatrixBlock<Self, Row, Col>(std::forward<Self>(self), 0, to, 0, to);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto ContinuousMatrix<Derived>::topRightCorner(this auto&& self, size_t toRow, size_t fromCol) noexcept {
        using Self = decltype(self);
        return ContinuousMatrixBlock<Self, Row, Col>(std::forward<Self>(self), 0, toRow, fromCol, self.getRow() - fromCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto ContinuousMatrix<Derived>::bottomLeftCorner(this auto&& self, size_t fromRow, size_t toCol) noexcept {
        using Self = decltype(self);
        return ContinuousMatrixBlock<Self, Row, Col>(std::forward<Self>(self), fromRow, self.getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto ContinuousMatrix<Derived>::bottomRightCorner(this auto&& self, size_t fromRow, size_t fromCol) noexcept {
        using Self = decltype(self);
        return ContinuousMatrixBlock<Self, Row, Col>(std::forward<Self>(self), fromRow, self.getRow() - fromRow, fromCol, self.getCol() - fromCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto ContinuousMatrix<Derived>::bottomRightCorner(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        return ContinuousMatrixBlock<Self, Row, Col>(std::forward<Self>(self), from, self.getRow() - from, from, self.getCol() - from);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto ContinuousMatrix<Derived>::block(this auto&& self, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept {
        using Self = decltype(self);
        return ContinuousMatrixBlock<Self, Row, Col>(std::forward<Self>(self), fromRow, rowCount, fromCol, colCount);
    }

    template<class Derived>
    auto ContinuousMatrix<Derived>::flatten() noexcept {
        return FlattenC<Derived>(Base::getDerived());
    }

    template<class Derived>
    const auto ContinuousMatrix<Derived>::flatten() const noexcept {
        return FlattenC<Derived>(Base::getConstCastDerived());
    }

    template<class Derived>
    void ContinuousMatrix<Derived>::zeros() {
        Base::getDerived().flatten().zeros();
    }

    template<class Derived>
    void ContinuousMatrix<Derived>::read(const T* __restrict p) {
        Base::getDerived().flatten().read(p);
    }

    template<class Derived>
    template<RNG R>
    void ContinuousMatrix<Derived>::random_uniform() {
        Base::getDerived().flatten().template random_uniform<R>();
    }

    template<class Derived>
    template<RNG R>
    void ContinuousMatrix<Derived>::random_normal() {
        Base::getDerived().flatten().template random_normal<R>();
    }

#ifdef PHYSICA_HDF5
    template<class Derived>
    const H5DataSet<2> ContinuousMatrix<Derived>::read(const H5Loc& loc, const char* name) {
        const auto dataset = loc.openDataSet<2>(name);
        const size_t maxMajor = dataset.getSize(0);
        const size_t maxMinor = dataset.getSize(1);
        Base::resize(Base::rowFromMajorMinor(maxMajor, maxMinor), Base::colFromMajorMinor(maxMajor, maxMinor));

        const auto memSpace = H5DataSpace<1>({maxMinor});
        auto fileSpace = H5DataSpace<2>({Base::getMaxMajor(), Base::getMaxMinor()});
        for (size_t major = 0; major < maxMajor; ++major) {
            fileSpace.selectHyperslab(H5S_SELECT_SET, {1, maxMinor}, {major, 0});
            if constexpr (isColMatrix)
                dataset.read(Base::getDerived().col(major).data(), ScalarType::dtype_hdf5(), memSpace, fileSpace);
            else
                dataset.read(Base::getDerived().row(major).data(), ScalarType::dtype_hdf5(), memSpace, fileSpace);
        }
        return dataset;
    }

    template<class Derived>
    H5DataSet<2> ContinuousMatrix<Derived>::write(H5Loc& loc, const char* name) const {
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
            if constexpr (isColMatrix)
                dataset.write(Base::getDerived().col(major).data(), ScalarType::dtype_hdf5(), memSpace, fileSpace);
            else
                dataset.write(Base::getDerived().row(major).data(), ScalarType::dtype_hdf5(), memSpace, fileSpace);
        }
        return std::cref(dataset);
    }
#endif

    template<class Derived>
    auto ContinuousMatrix<Derived>::data() noexcept {
        return Base::getDerived().data();
    }

    template<class Derived>
    auto ContinuousMatrix<Derived>::data() const noexcept {
        return Base::getDerived().data();
    }

    template<class Derived>
    auto ContinuousMatrix<Derived>::data_ptr(this auto&& self, size_t r, size_t c) noexcept {
        assert(r < self.getRow());
        assert(c < self.getCol());
        if constexpr (isRowMatrix)
            return self.data() + r * self.getCol() + c;
        else
            return self.data() + c * self.getRow() + r;
    }
}
