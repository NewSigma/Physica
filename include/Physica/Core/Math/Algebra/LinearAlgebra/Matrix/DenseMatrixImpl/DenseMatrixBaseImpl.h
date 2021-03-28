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

#include <cassert>

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
    typename DenseMatrixBase<Derived>::ScalarType DenseMatrixBase<Derived>::determinate() const {
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
    void DenseMatrixBase<Derived>::rowReduce(size_t r1, size_t r2, size_t elementIndex) {
        Derived& matrix = this->getDerived();
        assert(!matrix(r1, elementIndex).isZero());
        const size_t column = matrix.getColumn();
        const ScalarType factor = matrix(r2, elementIndex) / matrix(r1, elementIndex);
        for (size_t i = 0; i < column; ++i)
            matrix(r1, i) -= matrix(r2, i) * factor;
        assert(matrix(r2, elementIndex).isZero());
    }
    /* Operators */
    template<class Derived>
    Derived operator+(const DenseMatrixBase<Derived>& m1, const DenseMatrixBase<Derived>& m2) {
        Q_ASSERT(m1.getRow() == m2.getRow() && m1.getColumn() == m2.getColumn());
        Derived result(m1.getRow(), m1.getColumn());
        auto ite1 = m1.cbegin();
        auto ite2 = m2.cbegin();
        for(auto ite = result.begin(); ite != result.end(); ++ite, ++ite1, ++ite2)
            *ite = *ite1 + *ite2;
        return result;
    }

    template<class Derived>
    Derived operator-(const DenseMatrixBase<Derived>& m1, const DenseMatrixBase<Derived>& m2) {
        Q_ASSERT(m1.getRow() == m2.getRow() && m1.getColumn() == m2.getColumn());
        Derived result(m1.getRow(), m1.getColumn());
        auto ite1 = m1.cbegin();
        auto ite2 = m2.cbegin();
        for(auto ite = result.begin(); ite != result.end(); ++ite, ++ite1, ++ite2)
            *ite = *ite1 - *ite2;
        return result;
    }

    template<class Derived>
    Derived operator*(const DenseMatrixBase<Derived>& m, const MultiScalar& n) {
        Derived result(m.getRow(), m.getColumn());
        auto ite1 = m.cbegin();
        for(auto ite = result.begin(); ite != result.end(); ++ite, ++ite1)
            *ite = *ite1 * n;
        return result;
    }

    template<class Derived>
    Derived operator-(const DenseMatrixBase<Derived>& m) {
        Derived result(m.getRow(), m.getColumn());
        auto ite1 = m.cbegin();
        for(auto ite = result.begin(); ite != result.end(); ++ite, ++ite1)
            *ite = -(*ite1);
        return result;
    }
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class Derived>
    Derived reciprocal(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(reciprocal(m[i]), i);
        result.setLength(length);
        return result;
    }

    template<class Derived>
    Derived sqrt(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(sqrt(m[i]), i);
        result.setLength(length);
        return result;
    }

    template<class Derived>
    Derived factorial(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.allocate(factorial(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived ln(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(ln(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived log(const DenseMatrixBase<Derived>& m, const MultiScalar& a) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(log(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived exp(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(exp(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived pow(const DenseMatrixBase<Derived>& m, const MultiScalar& a) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(pow(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived cos(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(cos(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived sin(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(sin(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived tan(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(tan(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived sec(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(sec(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived csc(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(csc(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived cot(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(cot(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arccos(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arccos(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arcsin(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arcsin(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arctan(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arctan(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arcsec(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arcsec(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arccsc(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arccsc(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arccot(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arccot(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived cosh(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(cosh(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived sinh(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(sinh(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived tanh(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(tanh(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived sech(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(sech(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived csch(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(csch(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived coth(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(coth(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arccosh(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arccosh(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arcsinh(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arcsinh(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arctanh(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arctanh(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arcsech(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arcsech(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arccsch(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arccsch(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arccoth(const DenseMatrixBase<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arccoth(m[i]), i);
        result.setLength(length);
        return result;
    }
}