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
    inline ContinuousMatrix<Derived>& ContinuousMatrix<Derived>::operator=(const ContinuousMatrix<Derived>& obj) {
        Base::operator=(obj);
        return *this;
    }

    template<class Derived>
    inline ContinuousMatrix<Derived>& ContinuousMatrix<Derived>::operator=(ContinuousMatrix<Derived>&& obj) noexcept {
        Base::operator=(std::forward<Base>(obj));
        return *this;
    }

    template<class Derived>
    Derived& ContinuousMatrix<Derived>::operator=(const ScalarType& s) {
        const size_t maxMajor = Base::getMaxMajor();
        const size_t maxMinor = Base::getMaxMinor();
        if constexpr (isReverseDiff) {
            using TracerType = typename ScalarType::TracerType;
            TracerType::getInstance().reserve(maxMajor * maxMinor);
            for (size_t i = 0; i < maxMajor; ++i)
                for (size_t j = 0; j < maxMinor; ++j)
                    Base::refFromMajorMinor(i, j) = s.copy();
        }
        else {
            for (size_t i = 0; i < maxMajor; ++i)
                for (size_t j = 0; j < maxMinor; ++j)
                    Base::refFromMajorMinor(i, j) = s;
        }
        return Base::getDerived();
    }

    template<class Derived>
    inline typename ContinuousMatrix<Derived>::RowVector ContinuousMatrix<Derived>::row(size_t r) {
        const bool useSpecialization = ContinuousMatrix<Derived>::ColumnAtCompile == 1;
        if constexpr (useSpecialization)
            return {Base::getDerived(), r, 1, 0};
        else
            return {Base::getDerived(), r, 0, Base::getColumn()};
    }

    template<class Derived>
    inline const typename ContinuousMatrix<Derived>::RowVector ContinuousMatrix<Derived>::row(size_t r) const {
        const bool useSpecialization = ContinuousMatrix<Derived>::ColumnAtCompile == 1;
        if constexpr (useSpecialization)
            return {Base::getConstCastDerived(), r, 1, 0};
        else {
            using ResultType = ContinuousMatrixBlock<Derived, 1, ContinuousMatrix<Derived>::ColumnAtCompile>;
            ResultType result{Base::getConstCastDerived(), r, 0, Base::getColumn()};
            if constexpr (isRowMatrix && isReverseDiff)
                assert(result.checkContinuous() && "[Error]: Const matrix must be continuous because we cannot modify it and make it continuous");
            return result;
        }
    }

    template<class Derived>
    inline typename ContinuousMatrix<Derived>::ColVector ContinuousMatrix<Derived>::col(size_t c) {
        return {Base::getDerived(), 0, Base::getRow(), c};
    }

    template<class Derived>
    inline const typename ContinuousMatrix<Derived>::ColVector ContinuousMatrix<Derived>::col(size_t c) const {
        using ResultType = ContinuousMatrixBlock<Derived, ContinuousMatrix<Derived>::RowAtCompile, 1>;
        ResultType result{Base::getConstCastDerived(), 0, Base::getRow(), c};
        if constexpr (isColumnMatrix && isReverseDiff)
            assert(result.checkContinuous() && "[Error]: Const matrix must be continuous because we cannot modify it and make it continuous");
        return result;
    }

    template<class Derived>
    template<size_t Row>
    inline typename ContinuousMatrix<Derived>::template RowBlock<Row>
    ContinuousMatrix<Derived>::rows(size_t fromRow, size_t rowCount) {
        return {Base::getDerived(), fromRow, rowCount, 0, Base::getColumn()};
    }

    template<class Derived>
    template<size_t Row>
    inline const typename ContinuousMatrix<Derived>::template RowBlock<Row>
    ContinuousMatrix<Derived>::rows(size_t fromRow, size_t rowCount) const {
        return {Base::getConstCastDerived(), fromRow, rowCount, 0, Base::getColumn()};
    }

    template<class Derived>
    template<size_t Row>
    inline typename ContinuousMatrix<Derived>::template RowBlock<Row>
    ContinuousMatrix<Derived>::topRows(size_t to) {
        return {Base::getDerived(), 0, to, 0, Base::getColumn()};
    }

    template<class Derived>
    template<size_t Row>
    inline const typename ContinuousMatrix<Derived>::template RowBlock<Row>
    ContinuousMatrix<Derived>::topRows(size_t to) const {
        return {Base::getConstCastDerived(), 0, to, 0, Base::getColumn()};
    }

    template<class Derived>
    template<size_t Row>
    inline typename ContinuousMatrix<Derived>::template RowBlock<Row>
    ContinuousMatrix<Derived>::bottomRows(size_t from) {
        return {Base::getDerived(), from, Base::getRow() - from, 0, Base::getColumn()};
    }

    template<class Derived>
    template<size_t Row>
    inline const typename ContinuousMatrix<Derived>::template RowBlock<Row>
    ContinuousMatrix<Derived>::bottomRows(size_t from) const {
        return {Base::getConstCastDerived(), from, Base::getRow() - from, 0, Base::getColumn()};
    }

    template<class Derived>
    template<size_t Column>
    inline typename ContinuousMatrix<Derived>::template ColBlock<Column>
    ContinuousMatrix<Derived>::cols(size_t fromCol, size_t colCount) {
        return {Base::getDerived(), 0, Base::getRow(), fromCol, colCount};
    }

    template<class Derived>
    template<size_t Column>
    inline const typename ContinuousMatrix<Derived>::template ColBlock<Column>
    ContinuousMatrix<Derived>::cols(size_t fromCol, size_t colCount) const {
        return {Base::getConstCastDerived(), 0, Base::getRow(), fromCol, colCount};
    }

    template<class Derived>
    template<size_t Column>
    inline typename ContinuousMatrix<Derived>::template ColBlock<Column>
    ContinuousMatrix<Derived>::leftCols(size_t to) {
        return {Base::getDerived(), 0, Base::getRow(), 0, to};
    }

    template<class Derived>
    template<size_t Column>
    inline const typename ContinuousMatrix<Derived>::template ColBlock<Column>
    ContinuousMatrix<Derived>::leftCols(size_t to) const {
        return {Base::getConstCastDerived(), 0, Base::getRow(), 0, to};
    }

    template<class Derived>
    template<size_t Column>
    inline typename ContinuousMatrix<Derived>::template ColBlock<Column>
    ContinuousMatrix<Derived>::rightCols(size_t from) {
        return {Base::getDerived(), 0, Base::getRow(), from, Base::getColumn() - from};
    }

    template<class Derived>
    template<size_t Column>
    inline const typename ContinuousMatrix<Derived>::template ColBlock<Column>
    ContinuousMatrix<Derived>::rightCols(size_t from) const {
        return {Base::getConstCastDerived(), 0, Base::getRow(), from, Base::getColumn() - from};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    inline ContinuousMatrixBlock<Derived, Row, Column>
    ContinuousMatrix<Derived>::topLeftCorner(size_t toRow, size_t toCol) {
        return {Base::getDerived(), 0, toRow, 0, toCol};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    inline const ContinuousMatrixBlock<Derived, Row, Column>
    ContinuousMatrix<Derived>::topLeftCorner(size_t toRow, size_t toCol) const {
        return {Base::getConstCastDerived(), 0, toRow, 0, toCol};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    inline ContinuousMatrixBlock<Derived, Row, Column>
    ContinuousMatrix<Derived>::topLeftCorner(size_t to) {
        return {Base::getDerived(), 0, to, 0, to};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    inline const ContinuousMatrixBlock<Derived, Row, Column>
    ContinuousMatrix<Derived>::topLeftCorner(size_t to) const {
        return {Base::getConstCastDerived(), 0, to, 0, to};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    inline ContinuousMatrixBlock<Derived, Row, Column>
    ContinuousMatrix<Derived>::topRightCorner(size_t toRow, size_t fromCol) {
        return {Base::getDerived(), 0, toRow, fromCol, Base::getRow() - fromCol};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    inline const ContinuousMatrixBlock<Derived, Row, Column>
    ContinuousMatrix<Derived>::topRightCorner(size_t toRow, size_t fromCol) const {
        return {Base::getConstCastDerived(), 0, toRow, fromCol, Base::getRow() - fromCol};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    inline ContinuousMatrixBlock<Derived, Row, Column>
    ContinuousMatrix<Derived>::bottomLeftCorner(size_t fromRow, size_t toCol) {
        return {Base::getDerived(), fromRow, Base::getRow() - fromRow, 0, toCol};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    inline const ContinuousMatrixBlock<Derived, Row, Column>
    ContinuousMatrix<Derived>::bottomLeftCorner(size_t fromRow, size_t toCol) const {
        return {Base::getConstCastDerived(), fromRow, Base::getRow() - fromRow, 0, toCol};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    inline ContinuousMatrixBlock<Derived, Row, Column>
    ContinuousMatrix<Derived>::bottomRightCorner(size_t fromRow, size_t fromCol) {
        return {Base::getDerived(), fromRow, Base::getRow() - fromRow, fromCol, Base::getColumn() - fromCol};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    inline const ContinuousMatrixBlock<Derived, Row, Column>
    ContinuousMatrix<Derived>::bottomRightCorner(size_t fromRow, size_t fromCol) const {
        return {Base::getConstCastDerived(), fromRow, Base::getRow() - fromRow, fromCol, Base::getColumn() - fromCol};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    inline ContinuousMatrixBlock<Derived, Row, Column> ContinuousMatrix<Derived>::bottomRightCorner(size_t from) {
        return {Base::getDerived(), from, Base::getRow() - from, from, Base::getColumn() - from};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    inline const ContinuousMatrixBlock<Derived, Row, Column> ContinuousMatrix<Derived>::bottomRightCorner(size_t from) const {
        return {Base::getConstCastDerived(), from, Base::getRow() - from, from, Base::getColumn() - from};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    inline ContinuousMatrixBlock<Derived, Row, Column> ContinuousMatrix<Derived>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) {
        return {Base::getDerived(), fromRow, rowCount, fromCol, colCount};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    inline const ContinuousMatrixBlock<Derived, Row, Column> ContinuousMatrix<Derived>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) const {
        return {Base::getConstCastDerived(), fromRow, rowCount, fromCol, colCount};
    }

    template<class Derived>
    inline void ContinuousMatrix<Derived>::makeContinuous() {
        if constexpr (isReverseDiff)
            Base::getDerived() = Base::getDerived().copy();
    }

    template<class Derived>
    bool ContinuousMatrix<Derived>::checkContinuous() const {
        if constexpr (isReverseDiff) {
            if constexpr (MatrixOption::isRowMatrix<Derived>()) {
                for (size_t i = 0; i < Base::getRow(); ++i)
                    if (!row(i).checkContinuous())
                        return false;
            }
            else {
                for (size_t i = 0; i < Base::getColumn(); ++i)
                    if (!col(i).checkContinuous())
                        return false;
            }
        }
        return true;
    }

    template<class Derived>
    template<class RandomGenerator>
    void ContinuousMatrix<Derived>::random_uniform(RandomGenerator& gen) {
        if constexpr (isReverseDiff) {
            using TracerType = typename ScalarType::TracerType;
            const size_t maxMajor = Base::getMaxMajor();
            const size_t maxMinor = Base::getMaxMinor();
            TracerType::getInstance().reserve(maxMajor * maxMinor);
            for (size_t major = 0; major < maxMajor; ++major)
                for (size_t minor = 0; minor < maxMinor; ++minor)
                    Base::refFromMajorMinor(major, minor) = ScalarType::random_uniform(gen);
        }
        else
            Base::random_uniform(gen);
    }

    template<class Derived>
    template<class RandomGenerator>
    void ContinuousMatrix<Derived>::random_normal(RandomGenerator& gen) {
        if constexpr (isReverseDiff) {
            using TracerType = typename ScalarType::TracerType;
            const size_t maxMajor = Base::getMaxMajor();
            const size_t maxMinor = Base::getMaxMinor();
            TracerType::getInstance().reserve(maxMajor * maxMinor);
            for (size_t major = 0; major < maxMajor; ++major)
                for (size_t minor = 0; minor < maxMinor; ++minor)
                    Base::refFromMajorMinor(major, minor) = ScalarType::random_normal(gen);
        }
        else
            Base::random_normal(gen);
    }

    template<class Derived>
    template<class Distribution, class RandomGenerator>
    void ContinuousMatrix<Derived>::random_any(Distribution& dist, RandomGenerator& gen) {
        if constexpr (isReverseDiff) {
            using TracerType = typename ScalarType::TracerType;
            const size_t maxMajor = Base::getMaxMajor();
            const size_t maxMinor = Base::getMaxMinor();
            TracerType::getInstance().reserve(maxMajor * maxMinor);
            for (size_t major = 0; major < maxMajor; ++major)
                for (size_t minor = 0; minor < maxMinor; ++minor)
                    Base::refFromMajorMinor(major, minor) = ScalarType::random_any(dist, gen);
        }
        else
            Base::random_any(dist, gen);
    }
#ifdef PHYSICA_HDF5
    template<class Derived>
    const H5DataSet<2> ContinuousMatrix<Derived>::read(const H5Location& loc, const char* name) {
        const auto dataset = loc.openDataSet<2>(name);
        const size_t maxMajor = dataset.getSize(0);
        const size_t maxMinor = dataset.getSize(1);
        resize(Base::rowFromMajorMinor(maxMajor, maxMinor), Base::columnFromMajorMinor(maxMajor, maxMinor));

        const auto memSpace = H5DataSpace<1>({maxMinor});
        auto fileSpace = H5DataSpace<2>({Base::getMaxMajor(), Base::getMaxMinor()});
        for (size_t major = 0; major < maxMajor; ++major) {
            fileSpace.selectHyperslab(H5S_SELECT_SET, {1, maxMinor}, {major, 0});
            if constexpr (isColumnMatrix)
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
            if constexpr (isColumnMatrix)
                dataset.write(col(major).data(), ScalarType::getH5DataType(), memSpace, fileSpace);
            else
                dataset.write(row(major).data(), ScalarType::getH5DataType(), memSpace, fileSpace);
        }
        return std::cref(dataset);
    }
#endif
}
