/*
 * Copyright 2021-2022 WeiBo He.
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
    namespace Internal {
        /**
         * \tparam rank
         * The rank of matrix.
         */
        template<class Derived, size_t rank>
        class DeterminateImpl {
        public:
            static typename Derived::ScalarType run(const Derived& m) {
                //TODO
                assert(false);
            }
        };

        template<class Derived>
        class DeterminateImpl<Derived, 1> {
        public:
            static inline typename Derived::ScalarType run(const Derived& m) {
                return m(0, 0);
            }
        };

        template<class Derived>
        class DeterminateImpl<Derived, 2> {
        public:
            static inline typename Derived::ScalarType run(const Derived& m) {
                return m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0);
            }
        };

        template<class Derived>
        class DeterminateImpl<Derived, 3> {
        public:
            static inline typename Derived::ScalarType run(const Derived& m) {
            return m(0, 0) * (m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1))
                    + m(0, 1) * (m(1, 2) * m(2, 0) - m(1, 0) * m(2, 2))
                    + m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));
            }
        };
    }

    template<class Derived>
    template<class OtherMatrix>
    Derived& LValueMatrix<Derived>::operator=(const RValueMatrix<OtherMatrix>& m) {
        Base::getDerived().resize(m.getRow(), m.getColumn());
        m.assignTo(*this);
        return Base::getDerived();
    }
    
    template<class Derived>
    template<ScalarOption option, bool errorTrack>
    Derived& LValueMatrix<Derived>::operator=(const Scalar<option, errorTrack>& s) {
        for (size_t i = 0; i < getMaxMajor(); ++i)
            for (size_t j = 0; j < getMaxMinor(); ++j)
                getElementFromMajorMinor(i, j) = s;
        return Base::getDerived();
    }

    template<class Derived>
    inline typename LValueMatrix<Derived>::RowVector LValueMatrix<Derived>::row(size_t r) {
        return RowVector(Base::getDerived(), r, 0, Base::getColumn());
    }

    template<class Derived>
    inline const typename LValueMatrix<Derived>::RowVector LValueMatrix<Derived>::row(size_t r) const {
        return RowVector(Base::getConstCastDerived(), r, 0, Base::getColumn());
    }

    template<class Derived>
    inline typename LValueMatrix<Derived>::ColVector LValueMatrix<Derived>::col(size_t c) {
        return ColVector(Base::getDerived(), 0, Base::getRow(), c);
    }

    template<class Derived>
    inline const typename LValueMatrix<Derived>::ColVector LValueMatrix<Derived>::col(size_t c) const {
        return ColVector(Base::getConstCastDerived(), 0, Base::getRow(), c);
    }

    template<class Derived>
    inline MatrixBlock<Derived> LValueMatrix<Derived>::rows(size_t fromRow, size_t rowCount) {
        return MatrixBlock<Derived>(Base::getDerived(), fromRow, rowCount, 0, Base::getColumn());
    }

    template<class Derived>
    inline const MatrixBlock<Derived> LValueMatrix<Derived>::rows(size_t fromRow, size_t rowCount) const {
        return MatrixBlock<Derived>(Base::getConstCastDerived(), fromRow, rowCount, 0, Base::getColumn());
    }

    template<class Derived>
    inline MatrixBlock<Derived> LValueMatrix<Derived>::cols(size_t fromCol, size_t colCount) {
        return MatrixBlock<Derived>(Base::getDerived(), 0, Base::getRow(), fromCol, colCount);
    }

    template<class Derived>
    inline const MatrixBlock<Derived> LValueMatrix<Derived>::cols(size_t fromCol, size_t colCount) const {
        return MatrixBlock<Derived>(Base::getConstCastDerived(), 0, Base::getRow(), fromCol, colCount);
    }

    template<class Derived>
    inline MatrixBlock<Derived> LValueMatrix<Derived>::topRows(size_t to) {
        return MatrixBlock<Derived>(Base::getDerived(), 0, to, 0, Base::getColumn());
    }

    template<class Derived>
    inline const MatrixBlock<Derived> LValueMatrix<Derived>::topRows(size_t to) const {
        return MatrixBlock<Derived>(Base::getConstCastDerived(), 0, to, 0, Base::getColumn());
    }

    template<class Derived>
    inline MatrixBlock<Derived> LValueMatrix<Derived>::bottomRows(size_t from) {
        return MatrixBlock<Derived>(Base::getDerived(), from, Base::getRow() - from, 0, Base::getColumn());
    }

    template<class Derived>
    inline const MatrixBlock<Derived> LValueMatrix<Derived>::bottomRows(size_t from) const {
        return MatrixBlock<Derived>(Base::getConstCastDerived(), from, Base::getRow() - from, 0, Base::getColumn());
    }

    template<class Derived>
    inline MatrixBlock<Derived> LValueMatrix<Derived>::leftCols(size_t to) {
        return MatrixBlock<Derived>(Base::getDerived(), 0, Base::getRow(), 0, to);
    }

    template<class Derived>
    inline const MatrixBlock<Derived> LValueMatrix<Derived>::leftCols(size_t to) const {
        return MatrixBlock<Derived>(Base::getConstCastDerived(), 0, Base::getRow(), 0, to);
    }

    template<class Derived>
    inline MatrixBlock<Derived> LValueMatrix<Derived>::rightCols(size_t from) {
        return MatrixBlock<Derived>(Base::getDerived(), 0, Base::getRow(), from, Base::getColumn() - from);
    }

    template<class Derived>
    inline const MatrixBlock<Derived> LValueMatrix<Derived>::rightCols(size_t from) const {
        return MatrixBlock<Derived>(Base::getConstCastDerived(), 0, Base::getRow(), from, Base::getColumn() - from);
    }

    template<class Derived>
    inline MatrixBlock<Derived>
    LValueMatrix<Derived>::topLeftCorner(size_t toRow, size_t toCol) {
        return MatrixBlock<Derived>(Base::getDerived(), 0, toRow, 0, toCol);
    }

    template<class Derived>
    inline const MatrixBlock<Derived>
    LValueMatrix<Derived>::topLeftCorner(size_t toRow, size_t toCol) const {
        return MatrixBlock<Derived>(Base::getConstCastDerived(), 0, toRow, 0, toCol);
    }

    template<class Derived>
    inline MatrixBlock<Derived>
    LValueMatrix<Derived>::topLeftCorner(size_t to) {
        return MatrixBlock<Derived>(Base::getDerived(), 0, to, 0, to);
    }

    template<class Derived>
    inline const MatrixBlock<Derived>
    LValueMatrix<Derived>::topLeftCorner(size_t to) const {
        return MatrixBlock<Derived>(Base::getConstCastDerived(), 0, to, 0, to);
    }

    template<class Derived>
    inline MatrixBlock<Derived>
    LValueMatrix<Derived>::topRightCorner(size_t toRow, size_t fromCol) {
        return MatrixBlock<Derived>(Base::getDerived(), 0, toRow, fromCol, Base::getRow() - fromCol);
    }

    template<class Derived>
    inline const MatrixBlock<Derived>
    LValueMatrix<Derived>::topRightCorner(size_t toRow, size_t fromCol) const {
        return MatrixBlock<Derived>(Base::getConstCastDerived(), 0, toRow, fromCol, Base::getRow() - fromCol);
    }

    template<class Derived>
    inline MatrixBlock<Derived>
    LValueMatrix<Derived>::bottomLeftCorner(size_t fromRow, size_t toCol) {
        return MatrixBlock<Derived>(Base::getDerived(), fromRow, Base::getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    inline const MatrixBlock<Derived>
    LValueMatrix<Derived>::bottomLeftCorner(size_t fromRow, size_t toCol) const {
        return MatrixBlock<Derived>(Base::getConstCastDerived(), fromRow, Base::getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    inline MatrixBlock<Derived>
    LValueMatrix<Derived>::bottomRightCorner(size_t fromRow, size_t fromCol) {
        return MatrixBlock<Derived>(Base::getDerived(), fromRow, Base::getRow() - fromRow, fromCol, Base::getColumn() - fromCol);
    }

    template<class Derived>
    inline const MatrixBlock<Derived>
    LValueMatrix<Derived>::bottomRightCorner(size_t fromRow, size_t fromCol) const {
        return MatrixBlock<Derived>(Base::getConstCastDerived(), fromRow, Base::getRow() - fromRow, fromCol, Base::getColumn() - fromCol);
    }

    template<class Derived>
    inline MatrixBlock<Derived> LValueMatrix<Derived>::bottomRightCorner(size_t from) {
        return MatrixBlock<Derived>(Base::getDerived(), from, Base::getRow() - from, from, Base::getColumn() - from);
    }

    template<class Derived>
    inline const MatrixBlock<Derived> LValueMatrix<Derived>::bottomRightCorner(size_t from) const {
        return MatrixBlock<Derived>(Base::getConstCastDerived(), from, Base::getRow() - from, from, Base::getColumn() - from);
    }

    template<class Derived>
    inline MatrixBlock<Derived> LValueMatrix<Derived>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) {
        return MatrixBlock<Derived>(Base::getDerived(), fromRow, rowCount, fromCol, colCount);
    }

    template<class Derived>
    inline const MatrixBlock<Derived> LValueMatrix<Derived>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) const {
        return MatrixBlock<Derived>(Base::getConstCastDerived(), fromRow, rowCount, fromCol, colCount);
    }

    template<class Derived>
    InverseMatrix<Derived> LValueMatrix<Derived>::inverse() const noexcept {
        return InverseMatrix<Derived>(this->getDerived());
    }

    template<class Derived>
    typename LValueMatrix<Derived>::ScalarType LValueMatrix<Derived>::determinate() const {
        assert(Base::getDerived().getRow() == Base::getDerived().getColumn());
        using namespace Internal;
        constexpr size_t RowAtCompile = Traits<Derived>::RowAtCompile;
        constexpr size_t ColumnAtCompile = Traits<Derived>::ColumnAtCompile;
        //Row equals column at runtime from the assert, but RowAtCompile and ColumnAtCompile may not equal. Ether of them could be dynamic.
        constexpr size_t Rank = RowAtCompile > ColumnAtCompile ? RowAtCompile : ColumnAtCompile;
        return DeterminateImpl<Derived, Rank>::run(Base::getDerived());
    }
    /**
     * Reduce the element at one row using the other row.
     * \param r1
     * The index of row to be used.
     * \param r2
     * The index of row that the element belongs to.
     * \param elementIndex
     * Index of the element to be reduced.
     */
    template<class Derived>
    void LValueMatrix<Derived>::rowReduce(size_t r1, size_t r2, size_t elementIndex) {
        Derived& matrix = Base::getDerived();
        assert(!matrix(r1, elementIndex).isZero());
        const size_t column = matrix.getColumn();
        const ScalarType factor = matrix(r2, elementIndex) / matrix(r1, elementIndex);
        for (size_t i = 0; i < column; ++i)
            matrix(r1, i) -= matrix(r2, i) * factor;
        assert(matrix(r2, elementIndex).isZero());
    }

    template<class Derived>
    void LValueMatrix<Derived>::rowReduce(size_t r1, size_t r2, const ScalarType& factor) {
        Derived& matrix = Base::getDerived();
        const size_t column = matrix.getColumn();
        for (size_t i = 0; i < column; ++i)
            matrix(r1, i) -= matrix(r2, i) * factor;
    }

    template<class Derived>
    void LValueMatrix<Derived>::columnReduce(size_t c1, size_t c2, size_t elementIndex) {
        Derived& matrix = Base::getDerived();
        assert(!matrix(elementIndex, c1).isZero());
        const size_t row = matrix.getRow();
        const ScalarType factor = matrix(c2, elementIndex) / matrix(c1, elementIndex);
        for (size_t i = 0; i < row; ++i)
            matrix(i, c1) -= matrix(i, c2) * factor;
        assert(matrix(elementIndex, c2).isZero());
    }

    template<class Derived>
    void LValueMatrix<Derived>::columnReduce(size_t c1, size_t c2, const ScalarType& factor) {
        Derived& matrix = Base::getDerived();
        const size_t row = matrix.getRow();
        for (size_t i = 0; i < row; ++i)
            matrix(i, c1) -= matrix(i, c2) * factor;
    }

    template<class Derived>
    inline void LValueMatrix<Derived>::majorReduce(size_t v1, size_t v2, size_t elementIndex) {
        if constexpr (DenseMatrixOption::isColumnMatrix<Derived>())
            columnReduce(v1, v2, elementIndex);
        else
            rowReduce(v1, v2, elementIndex);
    }

    template<class Derived>
    inline void LValueMatrix<Derived>::majorReduce(size_t v1, size_t v2, const ScalarType& factor) {
        if constexpr (DenseMatrixOption::isColumnMatrix<Derived>())
            columnReduce(v1, v2, factor);
        else
            rowReduce(v1, v2, factor);
    }

    template<class Derived>
    void LValueMatrix<Derived>::rowMulScalar(size_t r, const ScalarType& factor) {
        const size_t column = Base::getColumn();
        for (size_t i = 0; i < column; ++i)
            (*this)(r, i) *= factor;
    }

    template<class Derived>
    void LValueMatrix<Derived>::columnMulScalar(size_t c, const ScalarType& factor) {
        const size_t row = Base::getRow();
        for (size_t i = 0; i < row; ++i)
            (*this)(i, c) *= factor;
    }

    template<class Derived>
    inline void LValueMatrix<Derived>::majorMulScalar(size_t v, const ScalarType& factor) {
        if constexpr (DenseMatrixOption::isColumnMatrix<Derived>())
            columnMulScalar(v, factor);
        else
            rowMulScalar(v, factor);
    }

    template<class Derived>
    inline void LValueMatrix<Derived>::majorSwap(size_t v1, size_t v2) {
        if constexpr (DenseMatrixOption::isColumnMatrix<Derived>())
            Base::getDerived().columnSwap(v1, v2);
        else
            Base::getDerived().rowSwap(v1, v2);
    }

    template<class Derived>
    template<class OtherDerived>
    void LValueMatrix<Derived>::assignTo(LValueMatrix<OtherDerived>& target) const {
        using TargetType = LValueMatrix<OtherDerived>;
        const size_t max_i = target.getMaxMajor();
        const size_t mat_j = target.getMaxMinor();
        for (size_t i = 0; i < max_i; ++i) {
            for (size_t j = 0; j < mat_j; ++j) {
                target.getElementFromMajorMinor(i, j) = (*this)(TargetType::rowFromMajorMinor(i, j),
                                                                TargetType::columnFromMajorMinor(i, j));
            }
        }
    }

    template<class Derived>
    typename LValueMatrix<Derived>::ScalarType LValueMatrix<Derived>::max() const {
        ScalarType result;
        if constexpr (isColumnMatrix) {
            result = col(0).max();
            for (size_t i = 1; i < Base::getColumn(); ++i) {
                ScalarType temp = col(i).max();
                if (temp > result)
                    result = temp;
            }
        }
        else {
            result = row(0).max();
            for (size_t i = 1; i < Base::getRow(); ++i) {
                ScalarType temp = row(i).max();
                if (temp > result)
                    result = temp;
            }
        }
        return result;
    }

    template<class Derived>
    typename LValueMatrix<Derived>::ScalarType LValueMatrix<Derived>::min() const {
        ScalarType result;
        if constexpr (isColumnMatrix) {
            result = col(0).max();
            for (size_t i = 1; i < Base::getColumn(); ++i) {
                ScalarType temp = col(i).max();
                if (temp < result)
                    result = temp;
            }
        }
        else {
            result = row(0).max();
            for (size_t i = 1; i < Base::getRow(); ++i) {
                ScalarType temp = row(i).max();
                if (temp < result)
                    result = temp;
            }
        }
        return result;
    }

    template<class Derived>
    typename LValueMatrix<Derived>::ScalarType LValueMatrix<Derived>::sum() const {
        ScalarType result;
        if constexpr (isColumnMatrix) {
            result = col(0).asVector().sum();
            for (size_t i = 1; i < Base::getColumn(); ++i)
                result += col(i).asVector().sum();
        }
        else {
            result = row(0).asVector().sum();
            for (size_t i = 1; i < Base::getRow(); ++i)
                result += row(i).asVector().sum();
        }
        return result;
    }

    template<class Derived>
    inline size_t LValueMatrix<Derived>::getMaxMajor() const noexcept {
        if constexpr (DenseMatrixOption::isColumnMatrix<Derived>())
            return Base::getColumn();
        else
            return Base::getRow();
    }

    template<class Derived>
    inline size_t LValueMatrix<Derived>::getMaxMinor() const noexcept {
        if constexpr (DenseMatrixOption::isColumnMatrix<Derived>())
            return Base::getRow();
        else
            return Base::getColumn();
    }

    template<class Derived>
    typename LValueMatrix<Derived>::ScalarType& LValueMatrix<Derived>::getElementFromMajorMinor(size_t major, size_t minor) {
        size_t r, c;
        if constexpr(DenseMatrixOption::isColumnMatrix<Derived>()) {
            c = major;
            r = minor;
        }
        else {
            r = major;
            c = minor;
        }
        assert(r < Base::getDerived().getRow() && c < Base::getDerived().getColumn());
        return Base::getDerived()(r, c);
    }

    template<class Derived>
    const typename LValueMatrix<Derived>::ScalarType& LValueMatrix<Derived>::getElementFromMajorMinor(size_t major, size_t minor) const {
        size_t r, c;
        if constexpr(DenseMatrixOption::isColumnMatrix<Derived>()) {
            c = major;
            r = minor;
        }
        else {
            r = major;
            c = minor;
        }
        assert(r < Base::getDerived().getRow() && c < Base::getDerived().getColumn());
        return Base::getDerived()(r, c);
    }

    template<class Derived>
    void LValueMatrix<Derived>::toUnitMatrix() {
        assert(Base::getRow() == Base::getColumn());
        const size_t order = Base::getRow();
        for (size_t i = 0; i < order; ++i)
            for (size_t j = 0; j < order; ++j)
                getElementFromMajorMinor(i, j) = i == j ? ScalarType(1) : ScalarType(0);
    }

    template<class Derived>
    inline size_t LValueMatrix<Derived>::rowFromMajorMinor([[maybe_unused]] size_t major, [[maybe_unused]] size_t minor) noexcept {
        if constexpr (DenseMatrixOption::isColumnMatrix<Derived>())
            return minor;
        else
            return major;
    }

    template<class Derived>
    inline size_t LValueMatrix<Derived>::columnFromMajorMinor([[maybe_unused]] size_t major, [[maybe_unused]] size_t minor) noexcept {
        if constexpr (DenseMatrixOption::isColumnMatrix<Derived>())
            return major;
        else
            return minor;
    }
}
