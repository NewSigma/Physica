/*
 * Copyright 2023 Weibo He.
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

namespace Physica::Core {
    template<class Derived>
    inline ContinuousMatrix<Derived>& ContinuousMatrix<Derived>::operator=(const This& obj) {
        Base::operator=(obj);
        return *this;
    }

    template<class Derived>
    Derived& ContinuousMatrix<Derived>::operator=(const ScalarType& s) {
        if constexpr (isReverseDiff) {
            using TracerType = ScalarType::TracerType;
            const size_t maxMajor = Base::getMaxMajor();
            const size_t maxMinor = Base::getMaxMinor();
            TracerType::getInstance().reserve(maxMajor * maxMinor);
        }
        return Base::operator=(s); // TODO: Optimize it using SIMD
    }

    template<class Derived>
    inline ContinuousMatrix<Derived>::RowVector ContinuousMatrix<Derived>::row(size_t r) {
        const bool useSpecialization = ContinuousMatrix<Derived>::ColAtCompile == 1;
        if constexpr (useSpecialization)
            return {Base::getDerived(), r, 1, 0};
        else
            return {Base::getDerived(), r, 0, Base::getCol()};
    }

    template<class Derived>
    inline const ContinuousMatrix<Derived>::RowVector ContinuousMatrix<Derived>::row(size_t r) const {
        const bool useSpecialization = ContinuousMatrix<Derived>::ColAtCompile == 1;
        if constexpr (useSpecialization)
            return {Base::getConstCastDerived(), r, 1, 0};
        else
            return RowVector(Base::getConstCastDerived(), r, 0, Base::getCol());
    }

    template<class Derived>
    inline ContinuousMatrix<Derived>::ColVector ContinuousMatrix<Derived>::col(size_t c) {
        return {Base::getDerived(), 0, Base::getRow(), c};
    }

    template<class Derived>
    inline const ContinuousMatrix<Derived>::ColVector ContinuousMatrix<Derived>::col(size_t c) const {
        return ColVector(Base::getConstCastDerived(), 0, Base::getRow(), c);
    }

    template<class Derived>
    template<size_t Row>
    inline ContinuousMatrix<Derived>::template RowBlock<Row>
    ContinuousMatrix<Derived>::rows(size_t fromRow, size_t rowCount) {
        return {Base::getDerived(), fromRow, rowCount, 0, Base::getCol()};
    }

    template<class Derived>
    template<size_t Row>
    inline const ContinuousMatrix<Derived>::template RowBlock<Row>
    ContinuousMatrix<Derived>::rows(size_t fromRow, size_t rowCount) const {
        return {Base::getConstCastDerived(), fromRow, rowCount, 0, Base::getCol()};
    }

    template<class Derived>
    template<size_t Row>
    inline ContinuousMatrix<Derived>::template RowBlock<Row>
    ContinuousMatrix<Derived>::topRows(size_t to) {
        return {Base::getDerived(), 0, to, 0, Base::getCol()};
    }

    template<class Derived>
    template<size_t Row>
    inline const ContinuousMatrix<Derived>::template RowBlock<Row>
    ContinuousMatrix<Derived>::topRows(size_t to) const {
        return {Base::getConstCastDerived(), 0, to, 0, Base::getCol()};
    }

    template<class Derived>
    template<size_t Row>
    inline ContinuousMatrix<Derived>::template RowBlock<Row>
    ContinuousMatrix<Derived>::bottomRows(size_t from) {
        return {Base::getDerived(), from, Base::getRow() - from, 0, Base::getCol()};
    }

    template<class Derived>
    template<size_t Row>
    inline const ContinuousMatrix<Derived>::template RowBlock<Row>
    ContinuousMatrix<Derived>::bottomRows(size_t from) const {
        return {Base::getConstCastDerived(), from, Base::getRow() - from, 0, Base::getCol()};
    }

    template<class Derived>
    template<size_t Col>
    inline ContinuousMatrix<Derived>::template ColBlock<Col>
    ContinuousMatrix<Derived>::cols(size_t fromCol, size_t colCount) {
        return {Base::getDerived(), 0, Base::getRow(), fromCol, colCount};
    }

    template<class Derived>
    template<size_t Col>
    inline const ContinuousMatrix<Derived>::template ColBlock<Col>
    ContinuousMatrix<Derived>::cols(size_t fromCol, size_t colCount) const {
        return {Base::getConstCastDerived(), 0, Base::getRow(), fromCol, colCount};
    }

    template<class Derived>
    template<size_t Col>
    inline ContinuousMatrix<Derived>::template ColBlock<Col>
    ContinuousMatrix<Derived>::leftCols(size_t to) {
        return {Base::getDerived(), 0, Base::getRow(), 0, to};
    }

    template<class Derived>
    template<size_t Col>
    inline const ContinuousMatrix<Derived>::template ColBlock<Col>
    ContinuousMatrix<Derived>::leftCols(size_t to) const {
        return {Base::getConstCastDerived(), 0, Base::getRow(), 0, to};
    }

    template<class Derived>
    template<size_t Col>
    inline ContinuousMatrix<Derived>::template ColBlock<Col>
    ContinuousMatrix<Derived>::rightCols(size_t from) {
        return {Base::getDerived(), 0, Base::getRow(), from, Base::getCol() - from};
    }

    template<class Derived>
    template<size_t Col>
    inline const ContinuousMatrix<Derived>::template ColBlock<Col>
    ContinuousMatrix<Derived>::rightCols(size_t from) const {
        return {Base::getConstCastDerived(), 0, Base::getRow(), from, Base::getCol() - from};
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    inline ContinuousMatrixBlock<Derived, Row, Col>
    ContinuousMatrix<Derived>::topLeftCorner(size_t toRow, size_t toCol) {
        return {Base::getDerived(), 0, toRow, 0, toCol};
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    inline const ContinuousMatrixBlock<Derived, Row, Col>
    ContinuousMatrix<Derived>::topLeftCorner(size_t toRow, size_t toCol) const {
        return {Base::getConstCastDerived(), 0, toRow, 0, toCol};
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    inline ContinuousMatrixBlock<Derived, Row, Col>
    ContinuousMatrix<Derived>::topLeftCorner(size_t to) {
        return {Base::getDerived(), 0, to, 0, to};
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    inline const ContinuousMatrixBlock<Derived, Row, Col>
    ContinuousMatrix<Derived>::topLeftCorner(size_t to) const {
        return {Base::getConstCastDerived(), 0, to, 0, to};
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    inline ContinuousMatrixBlock<Derived, Row, Col>
    ContinuousMatrix<Derived>::topRightCorner(size_t toRow, size_t fromCol) {
        return {Base::getDerived(), 0, toRow, fromCol, Base::getRow() - fromCol};
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    inline const ContinuousMatrixBlock<Derived, Row, Col>
    ContinuousMatrix<Derived>::topRightCorner(size_t toRow, size_t fromCol) const {
        return {Base::getConstCastDerived(), 0, toRow, fromCol, Base::getRow() - fromCol};
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    inline ContinuousMatrixBlock<Derived, Row, Col>
    ContinuousMatrix<Derived>::bottomLeftCorner(size_t fromRow, size_t toCol) {
        return {Base::getDerived(), fromRow, Base::getRow() - fromRow, 0, toCol};
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    inline const ContinuousMatrixBlock<Derived, Row, Col>
    ContinuousMatrix<Derived>::bottomLeftCorner(size_t fromRow, size_t toCol) const {
        return {Base::getConstCastDerived(), fromRow, Base::getRow() - fromRow, 0, toCol};
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    inline ContinuousMatrixBlock<Derived, Row, Col>
    ContinuousMatrix<Derived>::bottomRightCorner(size_t fromRow, size_t fromCol) {
        return {Base::getDerived(), fromRow, Base::getRow() - fromRow, fromCol, Base::getCol() - fromCol};
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    inline const ContinuousMatrixBlock<Derived, Row, Col>
    ContinuousMatrix<Derived>::bottomRightCorner(size_t fromRow, size_t fromCol) const {
        return {Base::getConstCastDerived(), fromRow, Base::getRow() - fromRow, fromCol, Base::getCol() - fromCol};
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    inline ContinuousMatrixBlock<Derived, Row, Col> ContinuousMatrix<Derived>::bottomRightCorner(size_t from) {
        return {Base::getDerived(), from, Base::getRow() - from, from, Base::getCol() - from};
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    inline const ContinuousMatrixBlock<Derived, Row, Col> ContinuousMatrix<Derived>::bottomRightCorner(size_t from) const {
        return {Base::getConstCastDerived(), from, Base::getRow() - from, from, Base::getCol() - from};
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    inline ContinuousMatrixBlock<Derived, Row, Col> ContinuousMatrix<Derived>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) {
        return {Base::getDerived(), fromRow, rowCount, fromCol, colCount};
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    inline const ContinuousMatrixBlock<Derived, Row, Col> ContinuousMatrix<Derived>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) const {
        return {Base::getConstCastDerived(), fromRow, rowCount, fromCol, colCount};
    }

    template<class Derived>
    template<RandomGenerator R>
    void ContinuousMatrix<Derived>::random_uniform() {
        if constexpr (isReverseDiff) {
            using TracerType = ScalarType::TracerType;
            const size_t maxMajor = Base::getMaxMajor();
            const size_t maxMinor = Base::getMaxMinor();
            TracerType::getInstance().reserve(maxMajor * maxMinor);
            for (size_t major = 0; major < maxMajor; ++major)
                for (size_t minor = 0; minor < maxMinor; ++minor)
                    Base::refFromMajorMinor(major, minor) = ScalarType::template random_uniform<R>();
        }
        else
            Base::template random_uniform<R>();
    }

    template<class Derived>
    template<RandomGenerator R>
    void ContinuousMatrix<Derived>::random_normal() {
        if constexpr (isReverseDiff) {
            using TracerType = ScalarType::TracerType;
            const size_t maxMajor = Base::getMaxMajor();
            const size_t maxMinor = Base::getMaxMinor();
            TracerType::getInstance().reserve(maxMajor * maxMinor);
            for (size_t major = 0; major < maxMajor; ++major)
                for (size_t minor = 0; minor < maxMinor; ++minor)
                    Base::refFromMajorMinor(major, minor) = ScalarType::template random_normal<R>();
        }
        else
            Base::template random_normal<R>();
    }

    template<class Derived>
    template<class Distribution, RandomGenerator R>
    void ContinuousMatrix<Derived>::random_any(Distribution& dist) {
        if constexpr (isReverseDiff) {
            using TracerType = ScalarType::TracerType;
            const size_t maxMajor = Base::getMaxMajor();
            const size_t maxMinor = Base::getMaxMinor();
            TracerType::getInstance().reserve(maxMajor * maxMinor);
            for (size_t major = 0; major < maxMajor; ++major)
                for (size_t minor = 0; minor < maxMinor; ++minor)
                    Base::refFromMajorMinor(major, minor) = ScalarType::template random_any<decltype(dist), R>(dist);
        }
        else
            Base::template random_any<R>(dist);
    }
#ifdef PHYSICA_HDF5
    template<class Derived>
    const H5DataSet<2> ContinuousMatrix<Derived>::read(const H5Location& loc, const char* name) {
        const auto dataset = loc.openDataSet<2>(name);
        const size_t maxMajor = dataset.getSize(0);
        const size_t maxMinor = dataset.getSize(1);
        resize(Base::rowFromMajorMinor(maxMajor, maxMinor), Base::colFromMajorMinor(maxMajor, maxMinor));

        const auto memSpace = H5DataSpace<1>({maxMinor});
        auto fileSpace = H5DataSpace<2>({Base::getMaxMajor(), Base::getMaxMinor()});
        for (size_t major = 0; major < maxMajor; ++major) {
            fileSpace.selectHyperslab(H5S_SELECT_SET, {1, maxMinor}, {major, 0});
            if constexpr (isColMatrix)
                dataset.read(col(major).data(), ScalarType::getH5DataType(), memSpace, fileSpace);
            else
                dataset.read(row(major).data(), ScalarType::getH5DataType(), memSpace, fileSpace);
        }
        return dataset;
    }

    template<class Derived>
    H5DataSet<2> ContinuousMatrix<Derived>::write(H5Location& loc, const char* name) const {
        const size_t maxMajor = Base::getMaxMajor();
        const size_t maxMinor = Base::getMaxMinor();
        const auto space = H5DataSpace<2>({maxMajor, maxMinor});
        H5DataSet<2> dataset;
        if (loc.exists(name))
            dataset = loc.openDataSet<2>(name);
        else
            dataset = loc.createDataSet<2>(name, ScalarType::getH5DataType(), space);

        const auto memSpace = H5DataSpace<1>({maxMinor});
        auto fileSpace = H5DataSpace<2>({Base::getMaxMajor(), Base::getMaxMinor()});
        for (size_t major = 0; major < maxMajor; ++major) {
            fileSpace.selectHyperslab(H5S_SELECT_SET, {1, maxMinor}, {major, 0});
            if constexpr (isColMatrix)
                dataset.write(col(major).data(), ScalarType::getH5DataType(), memSpace, fileSpace);
            else
                dataset.write(row(major).data(), ScalarType::getH5DataType(), memSpace, fileSpace);
        }
        return std::cref(dataset);
    }
#endif
}
