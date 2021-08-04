/*
 * Copyright 2021 WeiBo He.
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
            static typename Derived::ScalarType run(const Derived& m) {
                //TODO
                assert(false);
            }
        };

        template<class Derived>
        class DeterminateImpl<Derived, 1> {
            static inline typename Derived::ScalarType run(const Derived& m) {
                return m(0, 0);
            }
        };

        template<class Derived>
        class DeterminateImpl<Derived, 2> {
            static inline typename Derived::ScalarType run(const Derived& m) {
                return m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0);
            }
        };

        template<class Derived>
        class DeterminateImpl<Derived, 3> {
            static inline typename Derived::ScalarType run(const Derived& m) {
            return m(0, 0) * (m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1))
                    + m(0, 1) * (m(1, 2) * m(2, 0) - m(1, 0) * m(2, 2))
                    + m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));
            }
        };
    }

    template<class Derived>
    typename LValueMatrix<Derived>::RowVector LValueMatrix<Derived>::row(size_t r) {
        return RowVector(*this, r, 0, Base::getColumn());
    }

    template<class Derived>
    typename LValueMatrix<Derived>::ColVector LValueMatrix<Derived>::col(size_t c) {
        return ColVector(*this, c, 0, Base::getRow());
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
}
