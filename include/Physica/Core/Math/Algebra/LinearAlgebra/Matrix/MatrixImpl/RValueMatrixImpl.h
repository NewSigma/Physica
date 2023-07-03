/*
 * Copyright 2021-2023 WeiBo He.
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

#include <cassert>

namespace Physica::Core {
    template<class Derived>
    template<class OtherDerived>
    void RValueMatrix<Derived>::assignTo(LValueMatrix<OtherDerived>& target) const {
        assert(getRow() == target.getRow());
        assert(getColumn() == target.getColumn());
        for (size_t i = 0; i < target.getMaxMajor(); ++i)
            for (size_t j = 0; j < target.getMaxMinor(); ++j)
                target.refFromMajorMinor(i, j) = calc(target.rowFromMajorMinor(i, j), target.columnFromMajorMinor(i, j));
    }

    template<class Derived>
    inline typename RValueMatrix<Derived>::RowVector RValueMatrix<Derived>::row(size_t r) {
        return RowVector(Base::getDerived(), r, 0, getColumn());
    }

    template<class Derived>
    inline const typename RValueMatrix<Derived>::RowVector RValueMatrix<Derived>::row(size_t r) const {
        return RowVector(Base::getConstCastDerived(), r, 0, getColumn());
    }

    template<class Derived>
    inline typename RValueMatrix<Derived>::ColVector RValueMatrix<Derived>::col(size_t c) {
        return ColVector(Base::getDerived(), 0, getRow(), c);
    }

    template<class Derived>
    inline const typename RValueMatrix<Derived>::ColVector RValueMatrix<Derived>::col(size_t c) const {
        return ColVector(Base::getConstCastDerived(), 0, getRow(), c);
    }

    template<class Derived>
    inline RMatrixBlock<Derived> RValueMatrix<Derived>::rows(size_t fromRow, size_t rowCount) {
        return RMatrixBlock<Derived>(Base::getDerived(), fromRow, rowCount, 0, getColumn());
    }

    template<class Derived>
    inline const RMatrixBlock<Derived> RValueMatrix<Derived>::rows(size_t fromRow, size_t rowCount) const {
        return RMatrixBlock<Derived>(Base::getConstCastDerived(), fromRow, rowCount, 0, getColumn());
    }

    template<class Derived>
    inline RMatrixBlock<Derived> RValueMatrix<Derived>::topRows(size_t to) {
        return RMatrixBlock<Derived>(Base::getDerived(), 0, to, 0, getColumn());
    }

    template<class Derived>
    inline const RMatrixBlock<Derived> RValueMatrix<Derived>::topRows(size_t to) const {
        return RMatrixBlock<Derived>(Base::getConstCastDerived(), 0, to, 0, getColumn());
    }

    template<class Derived>
    inline RMatrixBlock<Derived> RValueMatrix<Derived>::bottomRows(size_t from) {
        return RMatrixBlock<Derived>(Base::getDerived(), from, getRow() - from, 0, getColumn());
    }

    template<class Derived>
    inline const RMatrixBlock<Derived> RValueMatrix<Derived>::bottomRows(size_t from) const {
        return RMatrixBlock<Derived>(Base::getConstCastDerived(), from, getRow() - from, 0, getColumn());
    }

    template<class Derived>
    inline RMatrixBlock<Derived> RValueMatrix<Derived>::cols(size_t fromCol, size_t colCount) {
        return RMatrixBlock<Derived>(Base::getDerived(), 0, getRow(), fromCol, colCount);
    }

    template<class Derived>
    inline const RMatrixBlock<Derived> RValueMatrix<Derived>::cols(size_t fromCol, size_t colCount) const {
        return RMatrixBlock<Derived>(Base::getConstCastDerived(), 0, getRow(), fromCol, colCount);
    }

    template<class Derived>
    inline RMatrixBlock<Derived> RValueMatrix<Derived>::leftCols(size_t to) {
        return RMatrixBlock<Derived>(Base::getDerived(), 0, getRow(), 0, to);
    }

    template<class Derived>
    inline const RMatrixBlock<Derived> RValueMatrix<Derived>::leftCols(size_t to) const {
        return RMatrixBlock<Derived>(Base::getConstCastDerived(), 0, getRow(), 0, to);
    }

    template<class Derived>
    inline RMatrixBlock<Derived> RValueMatrix<Derived>::rightCols(size_t from) {
        return RMatrixBlock<Derived>(Base::getDerived(), 0, getRow(), from, getColumn() - from);
    }

    template<class Derived>
    inline const RMatrixBlock<Derived> RValueMatrix<Derived>::rightCols(size_t from) const {
        return RMatrixBlock<Derived>(Base::getConstCastDerived(), 0, getRow(), from, getColumn() - from);
    }

    template<class Derived>
    inline RMatrixBlock<Derived>
    RValueMatrix<Derived>::topLeftCorner(size_t toRow, size_t toCol) {
        return RMatrixBlock<Derived>(Base::getDerived(), 0, toRow, 0, toCol);
    }

    template<class Derived>
    inline const RMatrixBlock<Derived>
    RValueMatrix<Derived>::topLeftCorner(size_t toRow, size_t toCol) const {
        return RMatrixBlock<Derived>(Base::getConstCastDerived(), 0, toRow, 0, toCol);
    }

    template<class Derived>
    inline RMatrixBlock<Derived>
    RValueMatrix<Derived>::topLeftCorner(size_t to) {
        return RMatrixBlock<Derived>(Base::getDerived(), 0, to, 0, to);
    }

    template<class Derived>
    inline const RMatrixBlock<Derived>
    RValueMatrix<Derived>::topLeftCorner(size_t to) const {
        return RMatrixBlock<Derived>(Base::getConstCastDerived(), 0, to, 0, to);
    }

    template<class Derived>
    inline RMatrixBlock<Derived>
    RValueMatrix<Derived>::topRightCorner(size_t toRow, size_t fromCol) {
        return RMatrixBlock<Derived>(Base::getDerived(), 0, toRow, fromCol, getRow() - fromCol);
    }

    template<class Derived>
    inline const RMatrixBlock<Derived>
    RValueMatrix<Derived>::topRightCorner(size_t toRow, size_t fromCol) const {
        return RMatrixBlock<Derived>(Base::getConstCastDerived(), 0, toRow, fromCol, getRow() - fromCol);
    }

    template<class Derived>
    inline RMatrixBlock<Derived>
    RValueMatrix<Derived>::bottomLeftCorner(size_t fromRow, size_t toCol) {
        return RMatrixBlock<Derived>(Base::getDerived(), fromRow, getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    inline const RMatrixBlock<Derived>
    RValueMatrix<Derived>::bottomLeftCorner(size_t fromRow, size_t toCol) const {
        return RMatrixBlock<Derived>(Base::getConstCastDerived(), fromRow, getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    inline RMatrixBlock<Derived>
    RValueMatrix<Derived>::bottomRightCorner(size_t fromRow, size_t fromCol) {
        return RMatrixBlock<Derived>(Base::getDerived(), fromRow, getRow() - fromRow, fromCol, getColumn() - fromCol);
    }

    template<class Derived>
    inline const RMatrixBlock<Derived>
    RValueMatrix<Derived>::bottomRightCorner(size_t fromRow, size_t fromCol) const {
        return RMatrixBlock<Derived>(Base::getConstCastDerived(), fromRow, getRow() - fromRow, fromCol, getColumn() - fromCol);
    }

    template<class Derived>
    inline RMatrixBlock<Derived> RValueMatrix<Derived>::bottomRightCorner(size_t from) {
        return RMatrixBlock<Derived>(Base::getDerived(), from, getRow() - from, from, getColumn() - from);
    }

    template<class Derived>
    inline const RMatrixBlock<Derived> RValueMatrix<Derived>::bottomRightCorner(size_t from) const {
        return RMatrixBlock<Derived>(Base::getConstCastDerived(), from, getRow() - from, from, getColumn() - from);
    }

    template<class Derived>
    inline RMatrixBlock<Derived> RValueMatrix<Derived>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) {
        return RMatrixBlock<Derived>(Base::getDerived(), fromRow, rowCount, fromCol, colCount);
    }

    template<class Derived>
    inline const RMatrixBlock<Derived> RValueMatrix<Derived>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) const {
        return RMatrixBlock<Derived>(Base::getConstCastDerived(), fromRow, rowCount, fromCol, colCount);
    }

    template<class Derived>
    inline DiagVector<Derived, false> RValueMatrix<Derived>::diag() {
        return DiagVector<Derived, false>(Base::getDerived());
    }

    template<class Derived>
    inline const DiagVector<Derived, false> RValueMatrix<Derived>::diag() const {
        return DiagVector<Derived, false>(Base::getConstCastDerived());
    }

    template<class Derived>
    typename RValueMatrix<Derived>::ScalarType RValueMatrix<Derived>::calcFromMajorMinor(size_t major, size_t minor) const {
        return calc(MatrixOption::rowFromMajorMinor<Derived>(major, minor), MatrixOption::columnFromMajorMinor<Derived>(major, minor));
    }

    template<class Derived>
    typename RValueMatrix<Derived>::ScalarType RValueMatrix<Derived>::trace() const {
        assert(getRow() == getColumn());
        ScalarType result = ScalarType::Zero();
        for (size_t i = 0; i < getRow(); ++i)
            result += calc(i, i);
        return result;
    }

    template<class Derived>
    Transpose<Derived> RValueMatrix<Derived>::transpose() const noexcept {
        return Transpose<Derived>(Base::getDerived());
    }

    template<class Derived>
    Conjugate<Derived> RValueMatrix<Derived>::conjugate() const noexcept {
        return Conjugate<Derived>(Base::getDerived());
    }

    template<class Derived>
    Flatten<Derived> RValueMatrix<Derived>::flatten() const noexcept {
        return Flatten<Derived>(Base::getDerived());
    }

    template<class MatrixType, class MatrixType2>
    bool matrixNear(const RValueMatrix<MatrixType>& m1, const RValueMatrix<MatrixType2>& m2, double precision) {
        assert(m1.getRow() == m2.getRow());
        assert(m1.getColumn() == m2.getColumn());
        for (size_t i = 0; i < m1.getColumn(); ++i)
            for (size_t j = 0; j < m1.getRow(); ++j)
                if (!scalarNear(m1.calc(j, i), m2.calc(j, i), precision))
                    return false;
        return true;
    }

    template<class Derived>
    std::ostream& operator<<(std::ostream& os, const RValueMatrix<Derived>& m) {
        const auto row = m.getRow();
        const auto column = m.getColumn();
        //10 is the max precision of double.
        os << std::setprecision(10);
        for(size_t i = 0; i < row; ++i) {
            for(size_t j = 0; j < column; ++j)
                os << m.calc(i, j) << '\t';
            os << '\n';
        }
        //6 is the default precision.
        return os << std::setprecision(6);
    }
}