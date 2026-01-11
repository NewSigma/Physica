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
    auto ContinuousMatrix<Derived>::row(size_t r) noexcept {
        const bool IsMat1x1 = Base::ColAtCompile == 1;
        if constexpr (IsMat1x1)
            return ContinuousMatrixBlock<Derived, 1, 1>(Base::getDerived(), r, 0);
        else {
            if constexpr (isRowMatrix)
                return ContinuousMatrixBlock<Derived, 1, ColAtCompile>(Base::getDerived(), r, 0, Base::getCol());
            else
                return LMatrixBlock<Derived, 1, ColAtCompile>(Base::getDerived(), r, 0, Base::getCol());
        }
    }

    template<class Derived>
    const auto ContinuousMatrix<Derived>::row(size_t r) const noexcept {
        return Base::getConstCastDerived().row(r);
    }

    template<class Derived>
    auto ContinuousMatrix<Derived>::col(size_t c) noexcept {
        const bool IsMat1x1 = Base::RowAtCompile == 1;
        if constexpr (IsMat1x1)
            return ContinuousMatrixBlock<Derived, 1, 1>(Base::getDerived(), 0, c);
        else {
            if constexpr (isColMatrix)
                return ContinuousMatrixBlock<Derived, RowAtCompile, 1>(Base::getDerived(), 0, Base::getRow(), c);
            else
                return LMatrixBlock<Derived, RowAtCompile, 1>(Base::getDerived(), 0, Base::getRow(), c);
        }
    }

    template<class Derived>
    const auto ContinuousMatrix<Derived>::col(size_t c) const noexcept {
        return Base::getConstCastDerived().col(c);
    }

    template<class Derived>
    template<size_t Row>
    auto ContinuousMatrix<Derived>::rows(size_t fromRow, size_t rowCount) noexcept {
        return RowBlock<Row>(Base::getDerived(), fromRow, rowCount, 0, Base::getCol());
    }

    template<class Derived>
    template<size_t Row>
    const auto ContinuousMatrix<Derived>::rows(size_t fromRow, size_t rowCount) const noexcept {
        return RowBlock<Row>(Base::getConstCastDerived(), fromRow, rowCount, 0, Base::getCol());
    }

    template<class Derived>
    template<size_t Row>
    auto ContinuousMatrix<Derived>::topRows(size_t to) noexcept {
        return RowBlock<Row>(Base::getDerived(), 0, to, 0, Base::getCol());
    }

    template<class Derived>
    template<size_t Row>
    const auto ContinuousMatrix<Derived>::topRows(size_t to) const noexcept {
        return RowBlock<Row>(Base::getConstCastDerived(), 0, to, 0, Base::getCol());
    }

    template<class Derived>
    template<size_t Row>
    auto ContinuousMatrix<Derived>::bottomRows(size_t from) noexcept {
        return RowBlock<Row>(Base::getDerived(), from, Base::getRow() - from, 0, Base::getCol());
    }

    template<class Derived>
    template<size_t Row>
    const auto ContinuousMatrix<Derived>::bottomRows(size_t from) const noexcept {
        return RowBlock<Row>(Base::getConstCastDerived(), from, Base::getRow() - from, 0, Base::getCol());
    }

    template<class Derived>
    template<size_t Col>
    auto ContinuousMatrix<Derived>::cols(size_t fromCol, size_t colCount) noexcept {
        return ColBlock<Col>(Base::getDerived(), 0, Base::getRow(), fromCol, colCount);
    }

    template<class Derived>
    template<size_t Col>
    const auto ContinuousMatrix<Derived>::cols(size_t fromCol, size_t colCount) const noexcept {
        return ColBlock<Col>(Base::getConstCastDerived(), 0, Base::getRow(), fromCol, colCount);
    }

    template<class Derived>
    template<size_t Col>
    auto ContinuousMatrix<Derived>::leftCols(size_t to) noexcept {
        return ColBlock<Col>(Base::getDerived(), 0, Base::getRow(), 0, to);
    }

    template<class Derived>
    template<size_t Col>
    const auto ContinuousMatrix<Derived>::leftCols(size_t to) const noexcept {
        return ColBlock<Col>(Base::getConstCastDerived(), 0, Base::getRow(), 0, to);
    }

    template<class Derived>
    template<size_t Col>
    auto ContinuousMatrix<Derived>::rightCols(size_t from) noexcept {
        return ColBlock<Col>(Base::getDerived(), 0, Base::getRow(), from, Base::getCol() - from);
    }

    template<class Derived>
    template<size_t Col>
    const auto ContinuousMatrix<Derived>::rightCols(size_t from) const noexcept {
        return ColBlock<Col>(Base::getConstCastDerived(), 0, Base::getRow(), from, Base::getCol() - from);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto ContinuousMatrix<Derived>::topLeftCorner(size_t toRow, size_t toCol) noexcept {
        return BlockType<Row, Col>(Base::getDerived(), 0, toRow, 0, toCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    const auto ContinuousMatrix<Derived>::topLeftCorner(size_t toRow, size_t toCol) const noexcept {
        return BlockType<Row, Col>(Base::getConstCastDerived(), 0, toRow, 0, toCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto ContinuousMatrix<Derived>::topLeftCorner(size_t to) noexcept {
        return BlockType<Row, Col>(Base::getDerived(), 0, to, 0, to);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    const auto ContinuousMatrix<Derived>::topLeftCorner(size_t to) const noexcept {
        return BlockType<Row, Col>(Base::getConstCastDerived(), 0, to, 0, to);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto ContinuousMatrix<Derived>::topRightCorner(size_t toRow, size_t fromCol) noexcept {
        return BlockType<Row, Col>(Base::getDerived(), 0, toRow, fromCol, Base::getRow() - fromCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    const auto ContinuousMatrix<Derived>::topRightCorner(size_t toRow, size_t fromCol) const noexcept {
        return BlockType<Row, Col>(Base::getConstCastDerived(), 0, toRow, fromCol, Base::getRow() - fromCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto ContinuousMatrix<Derived>::bottomLeftCorner(size_t fromRow, size_t toCol) noexcept {
        return BlockType<Row, Col>(Base::getDerived(), fromRow, Base::getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    const auto ContinuousMatrix<Derived>::bottomLeftCorner(size_t fromRow, size_t toCol) const noexcept {
        return BlockType<Row, Col>(Base::getConstCastDerived(), fromRow, Base::getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto ContinuousMatrix<Derived>::bottomRightCorner(size_t fromRow, size_t fromCol) noexcept {
        return BlockType<Row, Col>(Base::getDerived(), fromRow, Base::getRow() - fromRow, fromCol, Base::getCol() - fromCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    const auto ContinuousMatrix<Derived>::bottomRightCorner(size_t fromRow, size_t fromCol) const noexcept {
        return BlockType<Row, Col>(Base::getConstCastDerived(), fromRow, Base::getRow() - fromRow, fromCol, Base::getCol() - fromCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto ContinuousMatrix<Derived>::bottomRightCorner(size_t from) noexcept {
        return BlockType<Row, Col>(Base::getDerived(), from, Base::getRow() - from, from, Base::getCol() - from);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    const auto ContinuousMatrix<Derived>::bottomRightCorner(size_t from) const noexcept {
        return BlockType<Row, Col>(Base::getConstCastDerived(), from, Base::getRow() - from, from, Base::getCol() - from);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto ContinuousMatrix<Derived>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept {
        return BlockType<Row, Col>(Base::getDerived(), fromRow, rowCount, fromCol, colCount);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    const auto ContinuousMatrix<Derived>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) const noexcept {
        return BlockType<Row, Col>(Base::getConstCastDerived(), fromRow, rowCount, fromCol, colCount);
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
                dataset.read(col(major).data(), ScalarType::dtype_hdf5(), memSpace, fileSpace);
            else
                dataset.read(row(major).data(), ScalarType::dtype_hdf5(), memSpace, fileSpace);
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
                dataset.write(col(major).data(), ScalarType::dtype_hdf5(), memSpace, fileSpace);
            else
                dataset.write(row(major).data(), ScalarType::dtype_hdf5(), memSpace, fileSpace);
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
