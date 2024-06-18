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

#include "LValueFlatten.h"

namespace Physica::Core {
    namespace Internal {
        /**
         * \tparam rank
         * The rank of matrix.
         */
        template<class Derived, size_t rank>
        class DeterminateImpl {
        public:
            static typename Derived::ScalarType run([[maybe_unused]] const Derived& m) {
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
        static_assert(RowAtCompile == Dynamic || OtherMatrix::RowAtCompile == Dynamic || RowAtCompile == OtherMatrix::RowAtCompile, "[Error]: Incompatible row number");
        static_assert(ColumnAtCompile == Dynamic || OtherMatrix::ColumnAtCompile == Dynamic || ColumnAtCompile == OtherMatrix::ColumnAtCompile, "[Error]: Incompatible column number");
        static_assert(!(!isComplex && OtherMatrix::isComplex), "[Error]: Cannot assign a complex matrix to real matrix");
        Base::getDerived().resize(m.getRow(), m.getColumn());
        m.getDerived().assignTo(Base::getDerived());
        return Base::getDerived();
    }
    
    template<class Derived>
    Derived& LValueMatrix<Derived>::operator=(const ScalarType& s) {
        for (size_t i = 0; i < Base::getMaxMajor(); ++i)
            for (size_t j = 0; j < Base::getMaxMinor(); ++j)
                refFromMajorMinor(i, j) = ScalarType(s.getDerived());
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
    inline LMatrixBlock<Derived> LValueMatrix<Derived>::rows(size_t fromRow, size_t rowCount) {
        return LMatrixBlock<Derived>(Base::getDerived(), fromRow, rowCount, 0, Base::getColumn());
    }

    template<class Derived>
    inline const LMatrixBlock<Derived> LValueMatrix<Derived>::rows(size_t fromRow, size_t rowCount) const {
        return LMatrixBlock<Derived>(Base::getConstCastDerived(), fromRow, rowCount, 0, Base::getColumn());
    }

    template<class Derived>
    inline LMatrixBlock<Derived> LValueMatrix<Derived>::topRows(size_t to) {
        return LMatrixBlock<Derived>(Base::getDerived(), 0, to, 0, Base::getColumn());
    }

    template<class Derived>
    inline const LMatrixBlock<Derived> LValueMatrix<Derived>::topRows(size_t to) const {
        return LMatrixBlock<Derived>(Base::getConstCastDerived(), 0, to, 0, Base::getColumn());
    }

    template<class Derived>
    inline LMatrixBlock<Derived> LValueMatrix<Derived>::bottomRows(size_t from) {
        return LMatrixBlock<Derived>(Base::getDerived(), from, Base::getRow() - from, 0, Base::getColumn());
    }

    template<class Derived>
    inline const LMatrixBlock<Derived> LValueMatrix<Derived>::bottomRows(size_t from) const {
        return LMatrixBlock<Derived>(Base::getConstCastDerived(), from, Base::getRow() - from, 0, Base::getColumn());
    }

    template<class Derived>
    inline LMatrixBlock<Derived> LValueMatrix<Derived>::cols(size_t fromCol, size_t colCount) {
        return LMatrixBlock<Derived>(Base::getDerived(), 0, Base::getRow(), fromCol, colCount);
    }

    template<class Derived>
    inline const LMatrixBlock<Derived> LValueMatrix<Derived>::cols(size_t fromCol, size_t colCount) const {
        return LMatrixBlock<Derived>(Base::getConstCastDerived(), 0, Base::getRow(), fromCol, colCount);
    }

    template<class Derived>
    inline LMatrixBlock<Derived> LValueMatrix<Derived>::leftCols(size_t to) {
        return LMatrixBlock<Derived>(Base::getDerived(), 0, Base::getRow(), 0, to);
    }

    template<class Derived>
    inline const LMatrixBlock<Derived> LValueMatrix<Derived>::leftCols(size_t to) const {
        return LMatrixBlock<Derived>(Base::getConstCastDerived(), 0, Base::getRow(), 0, to);
    }

    template<class Derived>
    inline LMatrixBlock<Derived> LValueMatrix<Derived>::rightCols(size_t from) {
        return LMatrixBlock<Derived>(Base::getDerived(), 0, Base::getRow(), from, Base::getColumn() - from);
    }

    template<class Derived>
    inline const LMatrixBlock<Derived> LValueMatrix<Derived>::rightCols(size_t from) const {
        return LMatrixBlock<Derived>(Base::getConstCastDerived(), 0, Base::getRow(), from, Base::getColumn() - from);
    }

    template<class Derived>
    inline LMatrixBlock<Derived>
    LValueMatrix<Derived>::topLeftCorner(size_t toRow, size_t toCol) {
        return LMatrixBlock<Derived>(Base::getDerived(), 0, toRow, 0, toCol);
    }

    template<class Derived>
    inline const LMatrixBlock<Derived>
    LValueMatrix<Derived>::topLeftCorner(size_t toRow, size_t toCol) const {
        return LMatrixBlock<Derived>(Base::getConstCastDerived(), 0, toRow, 0, toCol);
    }

    template<class Derived>
    inline LMatrixBlock<Derived>
    LValueMatrix<Derived>::topLeftCorner(size_t to) {
        return LMatrixBlock<Derived>(Base::getDerived(), 0, to, 0, to);
    }

    template<class Derived>
    inline const LMatrixBlock<Derived>
    LValueMatrix<Derived>::topLeftCorner(size_t to) const {
        return LMatrixBlock<Derived>(Base::getConstCastDerived(), 0, to, 0, to);
    }

    template<class Derived>
    inline LMatrixBlock<Derived>
    LValueMatrix<Derived>::topRightCorner(size_t toRow, size_t fromCol) {
        return LMatrixBlock<Derived>(Base::getDerived(), 0, toRow, fromCol, Base::getRow() - fromCol);
    }

    template<class Derived>
    inline const LMatrixBlock<Derived>
    LValueMatrix<Derived>::topRightCorner(size_t toRow, size_t fromCol) const {
        return LMatrixBlock<Derived>(Base::getConstCastDerived(), 0, toRow, fromCol, Base::getRow() - fromCol);
    }

    template<class Derived>
    inline LMatrixBlock<Derived>
    LValueMatrix<Derived>::bottomLeftCorner(size_t fromRow, size_t toCol) {
        return LMatrixBlock<Derived>(Base::getDerived(), fromRow, Base::getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    inline const LMatrixBlock<Derived>
    LValueMatrix<Derived>::bottomLeftCorner(size_t fromRow, size_t toCol) const {
        return LMatrixBlock<Derived>(Base::getConstCastDerived(), fromRow, Base::getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    inline LMatrixBlock<Derived>
    LValueMatrix<Derived>::bottomRightCorner(size_t fromRow, size_t fromCol) {
        return LMatrixBlock<Derived>(Base::getDerived(), fromRow, Base::getRow() - fromRow, fromCol, Base::getColumn() - fromCol);
    }

    template<class Derived>
    inline const LMatrixBlock<Derived>
    LValueMatrix<Derived>::bottomRightCorner(size_t fromRow, size_t fromCol) const {
        return LMatrixBlock<Derived>(Base::getConstCastDerived(), fromRow, Base::getRow() - fromRow, fromCol, Base::getColumn() - fromCol);
    }

    template<class Derived>
    inline LMatrixBlock<Derived> LValueMatrix<Derived>::bottomRightCorner(size_t from) {
        return LMatrixBlock<Derived>(Base::getDerived(), from, Base::getRow() - from, from, Base::getColumn() - from);
    }

    template<class Derived>
    inline const LMatrixBlock<Derived> LValueMatrix<Derived>::bottomRightCorner(size_t from) const {
        return LMatrixBlock<Derived>(Base::getConstCastDerived(), from, Base::getRow() - from, from, Base::getColumn() - from);
    }

    template<class Derived>
    inline LMatrixBlock<Derived> LValueMatrix<Derived>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) {
        return LMatrixBlock<Derived>(Base::getDerived(), fromRow, rowCount, fromCol, colCount);
    }

    template<class Derived>
    inline const LMatrixBlock<Derived> LValueMatrix<Derived>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) const {
        return LMatrixBlock<Derived>(Base::getConstCastDerived(), fromRow, rowCount, fromCol, colCount);
    }

    template<class Derived>
    inline DiagVector<Derived, true> LValueMatrix<Derived>::diag() {
        return DiagVector<Derived, true>(Base::getDerived());
    }

    template<class Derived>
    inline const DiagVector<Derived, true> LValueMatrix<Derived>::diag() const {
        return DiagVector<Derived, true>(Base::getConstCastDerived());
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
        assert(abs(matrix(r1, elementIndex)) > ScalarType(std::numeric_limits<ScalarType>::min()));
        const ScalarType factor = matrix(r2, elementIndex) / matrix(r1, elementIndex);
        auto row2 = row(r2);
        row2 -= row(r1).asVector() * factor;
        matrix(r2, elementIndex) = 0;
    }

    template<class Derived>
    void LValueMatrix<Derived>::columnReduce(size_t c1, size_t c2, size_t elementIndex) {
        Derived& matrix = Base::getDerived();
        assert(abs(matrix(elementIndex, c1)) > ScalarType(std::numeric_limits<ScalarType>::min()));
        const ScalarType factor = matrix(c2, elementIndex) / matrix(c1, elementIndex);
        auto col2 = col(c2);
        col2 -= col(c1).asVector() * factor;
        matrix(elementIndex, c2) = 0;
    }

    template<class Derived>
    inline void LValueMatrix<Derived>::majorReduce(size_t v1, size_t v2, size_t elementIndex) {
        if constexpr (MatrixOption::isColumnMatrix<Derived>())
            columnReduce(v1, v2, elementIndex);
        else
            rowReduce(v1, v2, elementIndex);
    }

    template<class Derived>
    inline void LValueMatrix<Derived>::majorReduce(size_t v1, size_t v2, const ScalarType& factor) {
        if constexpr (MatrixOption::isColumnMatrix<Derived>()) {
            auto col1 = col(v1);
            col1 -= col(v2).asVector() * factor;
        }
        else {
            auto row1 = row(v1);
            row1 -= row(v2).asVector() * factor;
        }
    }

    template<class Derived>
    inline void LValueMatrix<Derived>::majorMulScalar(size_t v, const ScalarType& factor) {
        if constexpr (MatrixOption::isColumnMatrix<Derived>()) {
            auto c = col(v);
            c.asVector() *= factor;
        }
        else {
            auto r = row(v);
            r.asVector() *= factor;
        }
    }

    template<class Derived>
    inline void LValueMatrix<Derived>::majorSwap(size_t v1, size_t v2) {
        if constexpr (MatrixOption::isColumnMatrix<Derived>())
            Base::getDerived().columnSwap(v1, v2);
        else
            Base::getDerived().rowSwap(v1, v2);
    }

    template<class Derived>
    template<class RandomGenerator>
    void LValueMatrix<Derived>::random_uniform(RandomGenerator& gen) {
        const size_t maxMajor = Base::getMaxMajor();
        const size_t maxMinor = Base::getMaxMinor();
        for (size_t major = 0; major < maxMajor; ++major)
            for (size_t minor = 0; minor < maxMinor; ++minor)
                refFromMajorMinor(major, minor) = ScalarType::random_uniform(gen);
    }

    template<class Derived>
    template<class RandomGenerator>
    void LValueMatrix<Derived>::random_normal(RandomGenerator& gen) {
        const size_t maxMajor = Base::getMaxMajor();
        const size_t maxMinor = Base::getMaxMinor();
        for (size_t major = 0; major < maxMajor; ++major)
            for (size_t minor = 0; minor < maxMinor; ++minor)
                refFromMajorMinor(major, minor) = ScalarType::random_normal(gen);
    }

    template<class Derived>
    template<class Distribution, class RandomGenerator>
    void LValueMatrix<Derived>::random_any(Distribution& dist, RandomGenerator& gen) {
        const size_t maxMajor = Base::getMaxMajor();
        const size_t maxMinor = Base::getMaxMinor();
        for (size_t major = 0; major < maxMajor; ++major)
            for (size_t minor = 0; minor < maxMinor; ++minor)
                refFromMajorMinor(major, minor) = ScalarType::random_any(dist, gen);
    }

    template<class Derived>
    inline typename LValueMatrix<Derived>::ScalarType& LValueMatrix<Derived>::refFromMajorMinor(size_t major, size_t minor) {
        const size_t r = MatrixOption::rowFromMajorMinor<Derived>(major, minor);
        const size_t c = MatrixOption::columnFromMajorMinor<Derived>(major, minor);
        assert(r < Base::getDerived().getRow() && c < Base::getDerived().getColumn());
        return Base::getDerived()(r, c);
    }

    template<class Derived>
    inline const typename LValueMatrix<Derived>::ScalarType& LValueMatrix<Derived>::refFromMajorMinor(size_t major, size_t minor) const {
        const size_t r = MatrixOption::rowFromMajorMinor<Derived>(major, minor);
        const size_t c = MatrixOption::columnFromMajorMinor<Derived>(major, minor);
        assert(r < Base::getDerived().getRow() && c < Base::getDerived().getColumn());
        return Base::getDerived()(r, c);
    }

    template<class Derived>
    void LValueMatrix<Derived>::toUnitMatrix() {
        assert(Base::getRow() == Base::getColumn());
        const size_t order = Base::getRow();
        for (size_t i = 0; i < order; ++i)
            for (size_t j = 0; j < order; ++j)
                refFromMajorMinor(i, j) = ScalarType(i == j ? 1 : 0);
    }

    template<class Derived>
    LValueFlatten<Derived> LValueMatrix<Derived>::flatten() {
        return LValueFlatten<Derived>(*this);
    }

    template<class Derived>
    const LValueFlatten<Derived> LValueMatrix<Derived>::flatten() const {
        return LValueFlatten<Derived>(*this);
    }
}
