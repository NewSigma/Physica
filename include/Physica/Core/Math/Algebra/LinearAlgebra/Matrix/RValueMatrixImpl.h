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
    template<class Derived>
    template<class OtherDerived>
    void RValueMatrix<Derived>::assignTo(LValueMatrix<OtherDerived>& target) const {
        assert(getRow() == target.getRow() && getColumn() == target.getColumn());
        for (size_t i = 0; i < target.getMaxMajor(); ++i)
            for (size_t j = 0; j < target.getMaxMinor(); ++j)
                target.getElementFromMajorMinor(i, j) = calc(target.rowFromMajorMinor(i, j), target.columnFromMajorMinor(i, j));
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
    inline RMatrixBlock<Derived> RValueMatrix<Derived>::cols(size_t fromCol, size_t colCount) {
        return RMatrixBlock<Derived>(Base::getDerived(), 0, getRow(), fromCol, colCount);
    }

    template<class Derived>
    inline const RMatrixBlock<Derived> RValueMatrix<Derived>::cols(size_t fromCol, size_t colCount) const {
        return RMatrixBlock<Derived>(Base::getConstCastDerived(), 0, getRow(), fromCol, colCount);
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
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class Derived>
    Derived reciprocal(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(reciprocal(m[i]), i);
        result.setLength(length);
        return result;
    }

    template<class Derived>
    Derived sqrt(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(sqrt(m[i]), i);
        result.setLength(length);
        return result;
    }

    template<class Derived>
    Derived factorial(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.allocate(factorial(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived ln(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(ln(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived log(const RValueMatrix<Derived>& m, const MultiScalar& a) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(log(m[i]), i);
        result.setLength(length);
        return result;
    }

    template<class Derived>
    Derived pow(const RValueMatrix<Derived>& m, const MultiScalar& a) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(pow(m[i]), i);
        result.setLength(length);
        return result;
    }

    template<class Derived>
    Derived tan(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(tan(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived sec(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(sec(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived csc(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(csc(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived cot(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(cot(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arccos(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arccos(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arcsin(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arcsin(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arctan(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arctan(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arcsec(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arcsec(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arccsc(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arccsc(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arccot(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arccot(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived cosh(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(cosh(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived sinh(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(sinh(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived tanh(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(tanh(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived sech(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(sech(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived csch(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(csch(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived coth(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(coth(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arccosh(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arccosh(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arcsinh(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arcsinh(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arctanh(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arctanh(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arcsech(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arcsech(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arccsch(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arccsch(m[i]), i);
        result.setLength(length);
        return result;
    }
    
    template<class Derived>
    Derived arccoth(const RValueMatrix<Derived>& m) {
        const auto length = m.getLength();
        Derived result();
        for(size_t i = 0; i < length; ++i)
            result.init(arccoth(m[i]), i);
        result.setLength(length);
        return result;
    }
}