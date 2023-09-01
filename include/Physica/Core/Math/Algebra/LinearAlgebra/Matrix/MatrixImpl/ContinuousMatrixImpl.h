/*
 * Copyright 2023 WeiBo He.
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
    inline ContinuousMatrixBlock<Derived, 1, ContinuousMatrix<Derived>::ColumnAtCompile> ContinuousMatrix<Derived>::row(size_t r) {
        const bool useSpecialization = ContinuousMatrix<Derived>::ColumnAtCompile == 1;
        if constexpr (useSpecialization)
            return {Base::getDerived(), r, 1, 0};
        else
            return {Base::getDerived(), r, 0, Base::getColumn()};
    }

    template<class Derived>
    inline const ContinuousMatrixBlock<Derived, 1, ContinuousMatrix<Derived>::ColumnAtCompile> ContinuousMatrix<Derived>::row(size_t r) const {
        const bool useSpecialization = ContinuousMatrix<Derived>::ColumnAtCompile == 1;
        if constexpr (useSpecialization)
            return {Base::getConstCastDerived(), r, 1, 0};
        else
            return {Base::getConstCastDerived(), r, 0, Base::getColumn()};
    }

    template<class Derived>
    inline ContinuousMatrixBlock<Derived, ContinuousMatrix<Derived>::RowAtCompile, 1> ContinuousMatrix<Derived>::col(size_t c) {
        return {Base::getDerived(), 0, Base::getRow(), c};
    }

    template<class Derived>
    inline const ContinuousMatrixBlock<Derived, ContinuousMatrix<Derived>::RowAtCompile, 1> ContinuousMatrix<Derived>::col(size_t c) const {
        return {Base::getConstCastDerived(), 0, Base::getRow(), c};
    }

    template<class Derived>
    template<size_t Row>
    inline ContinuousMatrixBlock<Derived, Row, ContinuousMatrix<Derived>::ColumnAtCompile> ContinuousMatrix<Derived>::rows(size_t fromRow, size_t rowCount) {
        return {Base::getDerived(), fromRow, rowCount, 0, Base::getColumn()};
    }

    template<class Derived>
    template<size_t Row>
    inline const ContinuousMatrixBlock<Derived, Row, ContinuousMatrix<Derived>::ColumnAtCompile> ContinuousMatrix<Derived>::rows(size_t fromRow, size_t rowCount) const {
        return {Base::getConstCastDerived(), fromRow, rowCount, 0, Base::getColumn()};
    }

    template<class Derived>
    template<size_t Row>
    inline ContinuousMatrixBlock<Derived, Row, ContinuousMatrix<Derived>::ColumnAtCompile>ContinuousMatrix<Derived>::topRows(size_t to) {
        return {Base::getDerived(), 0, to, 0, Base::getColumn()};
    }

    template<class Derived>
    template<size_t Row>
    inline const ContinuousMatrixBlock<Derived, Row, ContinuousMatrix<Derived>::ColumnAtCompile> ContinuousMatrix<Derived>::topRows(size_t to) const {
        return {Base::getConstCastDerived(), 0, to, 0, Base::getColumn()};
    }

    template<class Derived>
    template<size_t Row>
    inline ContinuousMatrixBlock<Derived, Row, ContinuousMatrix<Derived>::ColumnAtCompile>
    ContinuousMatrix<Derived>::bottomRows(size_t from) {
        return {Base::getDerived(), from, Base::getRow() - from, 0, Base::getColumn()};
    }

    template<class Derived>
    template<size_t Row>
    inline const ContinuousMatrixBlock<Derived, Row, ContinuousMatrix<Derived>::ColumnAtCompile>
    ContinuousMatrix<Derived>::bottomRows(size_t from) const {
        return {Base::getConstCastDerived(), from, Base::getRow() - from, 0, Base::getColumn()};
    }

    template<class Derived>
    template<size_t Column>
    inline ContinuousMatrixBlock<Derived, ContinuousMatrix<Derived>::RowAtCompile, Column> ContinuousMatrix<Derived>::cols(size_t fromCol, size_t colCount) {
        return {Base::getDerived(), 0, Base::getRow(), fromCol, colCount};
    }

    template<class Derived>
    template<size_t Column>
    inline const ContinuousMatrixBlock<Derived, ContinuousMatrix<Derived>::RowAtCompile, Column> ContinuousMatrix<Derived>::cols(size_t fromCol, size_t colCount) const {
        return {Base::getConstCastDerived(), 0, Base::getRow(), fromCol, colCount};
    }

    template<class Derived>
    template<size_t Column>
    inline ContinuousMatrixBlock<Derived, ContinuousMatrix<Derived>::RowAtCompile, Column> ContinuousMatrix<Derived>::leftCols(size_t to) {
        return {Base::getDerived(), 0, Base::getRow(), 0, to};
    }

    template<class Derived>
    template<size_t Column>
    inline const ContinuousMatrixBlock<Derived, ContinuousMatrix<Derived>::RowAtCompile, Column>
    ContinuousMatrix<Derived>::leftCols(size_t to) const {
        return {Base::getConstCastDerived(), 0, Base::getRow(), 0, to};
    }

    template<class Derived>
    template<size_t Column>
    inline ContinuousMatrixBlock<Derived, ContinuousMatrix<Derived>::RowAtCompile, Column>
    ContinuousMatrix<Derived>::rightCols(size_t from) {
        return {Base::getDerived(), 0, Base::getRow(), from, Base::getColumn() - from};
    }

    template<class Derived>
    template<size_t Column>
    inline const ContinuousMatrixBlock<Derived, ContinuousMatrix<Derived>::RowAtCompile, Column>
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
}
