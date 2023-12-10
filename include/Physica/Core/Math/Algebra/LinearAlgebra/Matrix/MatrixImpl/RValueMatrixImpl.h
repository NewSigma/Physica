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
        constexpr size_t OtherRow = OtherDerived::RowAtCompile;
        constexpr size_t OtherColumn = OtherDerived::ColumnAtCompile;
        static_assert(RowAtCompile == OtherRow || RowAtCompile == Dynamic || OtherRow == Dynamic, "Row mismatch between two matrix");
        static_assert(ColumnAtCompile == OtherColumn || ColumnAtCompile == Dynamic || OtherColumn == Dynamic, "Column mismatch between two matrix");
        using OtherScalar = typename OtherDerived::ScalarType;
        assert(getRow() == target.getRow());
        assert(getColumn() == target.getColumn());
        for (size_t i = 0; i < target.getMaxMajor(); ++i)
            for (size_t j = 0; j < target.getMaxMinor(); ++j)
                target.refFromMajorMinor(i, j) = OtherScalar(calc(target.rowFromMajorMinor(i, j), target.columnFromMajorMinor(i, j)));

        constexpr bool isContinuous = std::is_base_of<ContinuousMatrix<OtherDerived>, OtherDerived>::value;
        if constexpr (isContinuous && isReverseDiff)
            target.getDerived().makeContinuous();
    }

    template<class Derived>
    inline typename RValueMatrix<Derived>::RowVector RValueMatrix<Derived>::row(size_t r) {
        return {Base::getDerived(), r, 0, getColumn()};
    }

    template<class Derived>
    inline const typename RValueMatrix<Derived>::RowVector RValueMatrix<Derived>::row(size_t r) const {
        return {Base::getConstCastDerived(), r, 0, getColumn()};
    }

    template<class Derived>
    inline typename RValueMatrix<Derived>::ColVector RValueMatrix<Derived>::col(size_t c) {
        return {Base::getDerived(), 0, getRow(), c};
    }

    template<class Derived>
    inline const typename RValueMatrix<Derived>::ColVector RValueMatrix<Derived>::col(size_t c) const {
        return {Base::getConstCastDerived(), 0, getRow(), c};
    }

    template<class Derived>
    inline RMatrixBlock<Derived> RValueMatrix<Derived>::rows(size_t fromRow, size_t rowCount) {
        return {Base::getDerived(), fromRow, rowCount, 0, getColumn()};
    }

    template<class Derived>
    inline const RMatrixBlock<Derived> RValueMatrix<Derived>::rows(size_t fromRow, size_t rowCount) const {
        return {Base::getConstCastDerived(), fromRow, rowCount, 0, getColumn()};
    }

    template<class Derived>
    inline RMatrixBlock<Derived> RValueMatrix<Derived>::topRows(size_t to) {
        return {Base::getDerived(), 0, to, 0, getColumn()};
    }

    template<class Derived>
    inline const RMatrixBlock<Derived> RValueMatrix<Derived>::topRows(size_t to) const {
        return {Base::getConstCastDerived(), 0, to, 0, getColumn()};
    }

    template<class Derived>
    inline RMatrixBlock<Derived> RValueMatrix<Derived>::bottomRows(size_t from) {
        return {Base::getDerived(), from, getRow() - from, 0, getColumn()};
    }

    template<class Derived>
    inline const RMatrixBlock<Derived> RValueMatrix<Derived>::bottomRows(size_t from) const {
        return {Base::getConstCastDerived(), from, getRow() - from, 0, getColumn()};
    }

    template<class Derived>
    inline RMatrixBlock<Derived> RValueMatrix<Derived>::cols(size_t fromCol, size_t colCount) {
        return {Base::getDerived(), 0, getRow(), fromCol, colCount};
    }

    template<class Derived>
    inline const RMatrixBlock<Derived> RValueMatrix<Derived>::cols(size_t fromCol, size_t colCount) const {
        return {Base::getConstCastDerived(), 0, getRow(), fromCol, colCount};
    }

    template<class Derived>
    inline RMatrixBlock<Derived> RValueMatrix<Derived>::leftCols(size_t to) {
        return {Base::getDerived(), 0, getRow(), 0, to};
    }

    template<class Derived>
    inline const RMatrixBlock<Derived> RValueMatrix<Derived>::leftCols(size_t to) const {
        return {Base::getConstCastDerived(), 0, getRow(), 0, to};
    }

    template<class Derived>
    inline RMatrixBlock<Derived> RValueMatrix<Derived>::rightCols(size_t from) {
        return {Base::getDerived(), 0, getRow(), from, getColumn() - from};
    }

    template<class Derived>
    inline const RMatrixBlock<Derived> RValueMatrix<Derived>::rightCols(size_t from) const {
        return {Base::getConstCastDerived(), 0, getRow(), from, getColumn() - from};
    }

    template<class Derived>
    inline RMatrixBlock<Derived>
    RValueMatrix<Derived>::topLeftCorner(size_t toRow, size_t toCol) {
        return {Base::getDerived(), 0, toRow, 0, toCol};
    }

    template<class Derived>
    inline const RMatrixBlock<Derived>
    RValueMatrix<Derived>::topLeftCorner(size_t toRow, size_t toCol) const {
        return {Base::getConstCastDerived(), 0, toRow, 0, toCol};
    }

    template<class Derived>
    inline RMatrixBlock<Derived>
    RValueMatrix<Derived>::topLeftCorner(size_t to) {
        return {Base::getDerived(), 0, to, 0, to};
    }

    template<class Derived>
    inline const RMatrixBlock<Derived>
    RValueMatrix<Derived>::topLeftCorner(size_t to) const {
        return {Base::getConstCastDerived(), 0, to, 0, to};
    }

    template<class Derived>
    inline RMatrixBlock<Derived>
    RValueMatrix<Derived>::topRightCorner(size_t toRow, size_t fromCol) {
        return {Base::getDerived(), 0, toRow, fromCol, getRow() - fromCol};
    }

    template<class Derived>
    inline const RMatrixBlock<Derived>
    RValueMatrix<Derived>::topRightCorner(size_t toRow, size_t fromCol) const {
        return {Base::getConstCastDerived(), 0, toRow, fromCol, getRow() - fromCol};
    }

    template<class Derived>
    inline RMatrixBlock<Derived>
    RValueMatrix<Derived>::bottomLeftCorner(size_t fromRow, size_t toCol) {
        return {Base::getDerived(), fromRow, getRow() - fromRow, 0, toCol};
    }

    template<class Derived>
    inline const RMatrixBlock<Derived>
    RValueMatrix<Derived>::bottomLeftCorner(size_t fromRow, size_t toCol) const {
        return {Base::getConstCastDerived(), fromRow, getRow() - fromRow, 0, toCol};
    }

    template<class Derived>
    inline RMatrixBlock<Derived>
    RValueMatrix<Derived>::bottomRightCorner(size_t fromRow, size_t fromCol) {
        return {Base::getDerived(), fromRow, getRow() - fromRow, fromCol, getColumn() - fromCol};
    }

    template<class Derived>
    inline const RMatrixBlock<Derived>
    RValueMatrix<Derived>::bottomRightCorner(size_t fromRow, size_t fromCol) const {
        return {Base::getConstCastDerived(), fromRow, getRow() - fromRow, fromCol, getColumn() - fromCol};
    }

    template<class Derived>
    inline RMatrixBlock<Derived> RValueMatrix<Derived>::bottomRightCorner(size_t from) {
        return {Base::getDerived(), from, getRow() - from, from, getColumn() - from};
    }

    template<class Derived>
    inline const RMatrixBlock<Derived> RValueMatrix<Derived>::bottomRightCorner(size_t from) const {
        return {Base::getConstCastDerived(), from, getRow() - from, from, getColumn() - from};
    }

    template<class Derived>
    inline RMatrixBlock<Derived> RValueMatrix<Derived>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) {
        return {Base::getDerived(), fromRow, rowCount, fromCol, colCount};
    }

    template<class Derived>
    inline const RMatrixBlock<Derived> RValueMatrix<Derived>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) const {
        return {Base::getConstCastDerived(), fromRow, rowCount, fromCol, colCount};
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
    inline typename RValueMatrix<Derived>::ScalarType RValueMatrix<Derived>::calcFromMajorMinor(size_t major, size_t minor) const {
        return calc(MatrixOption::rowFromMajorMinor<Derived>(major, minor), MatrixOption::columnFromMajorMinor<Derived>(major, minor));
    }

    template<class Derived>
    typename RValueMatrix<Derived>::ScalarType RValueMatrix<Derived>::max() const {
        ScalarType result;
        if constexpr (isColumnMatrix) {
            result = col(0).max();
            for (size_t i = 1; i < getColumn(); ++i) {
                ScalarType temp = col(i).max();
                if (temp > result)
                    result = temp;
            }
        }
        else {
            result = row(0).max();
            for (size_t i = 1; i < getRow(); ++i) {
                ScalarType temp = row(i).max();
                if (temp > result)
                    result = temp;
            }
        }
        return result;
    }

    template<class Derived>
    typename RValueMatrix<Derived>::ScalarType RValueMatrix<Derived>::min() const {
        ScalarType result;
        if constexpr (isColumnMatrix) {
            result = col(0).min();
            for (size_t i = 1; i < getColumn(); ++i) {
                ScalarType temp = col(i).min();
                if (temp < result)
                    result = temp;
            }
        }
        else {
            result = row(0).min();
            for (size_t i = 1; i < getRow(); ++i) {
                ScalarType temp = row(i).min();
                if (temp < result)
                    result = temp;
            }
        }
        return result;
    }

    template<class Derived>
    typename RValueMatrix<Derived>::ScalarType RValueMatrix<Derived>::trace() const {
        assert(getRow() == getColumn());
        ScalarType result = ScalarType(0);
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
    RValueFlatten<Derived> RValueMatrix<Derived>::flatten() const noexcept {
        return RValueFlatten<Derived>(Base::getDerived());
    }

    template<class Derived>
    typename RValueMatrix<Derived>::ScalarType RValueMatrix<Derived>::sum() const {
        ScalarType result = 0;
        for (size_t major = 0; major < getMaxMajor(); ++major)
            for (size_t minor = 0; minor < getMaxMinor(); ++minor)
                result += calcFromMajorMinor(major, minor);
        return result;
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

    template<class Derived>
    bool operator==(const RValueMatrix<Derived>& m1, const RValueMatrix<Derived>& m2) {
        if (m1.getRow() != m2.getRow())
            return false;
        if (m1.getColumn() != m2.getColumn())
            return false;
        for (size_t major = 0; major < m1.getMaxMajor(); ++major)
            for (size_t minor = 0; minor < m1.getMaxMinor(); ++minor)
                if (m1.calcFromMajorMinor(major, minor) != m2.calcFromMajorMinor(major, minor))
                    return false;
        return true;
    }
}