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

#include "DenseMatrixImpl/DenseMatrixStorage/DenseMatrixStorage.h"
#include "RMatrixBlock.h"

namespace Physica::Core {
    template<class Derived> class LValueMatrix;
    template<class MatrixType> class Transpose;
    template<class MatrixType> class Conjugate;
    /**
     * The \class DenseRValueMatrix provide algorithms that a matrix should support.
     * 
     * \tparam Derived
     * A class that contains data structure for a matrix.
     */
    template<class Derived>
    class RValueMatrix : public Utils::CRTPBase<Derived> {
        using Base = Utils::CRTPBase<Derived>;
    public:
        using ScalarType = typename Internal::Traits<Derived>::ScalarType;
        constexpr static size_t RowAtCompile = Internal::Traits<Derived>::RowAtCompile;
        constexpr static size_t ColumnAtCompile = Internal::Traits<Derived>::ColumnAtCompile;
        constexpr static size_t MaxRowAtCompile = Internal::Traits<Derived>::MaxRowAtCompile;
        constexpr static size_t MaxColumnAtCompile = Internal::Traits<Derived>::MaxColumnAtCompile;
        constexpr static size_t SizeAtCompile = Internal::Traits<Derived>::SizeAtCompile;
        constexpr static size_t MaxSizeAtCompile = Internal::Traits<Derived>::MaxSizeAtCompile;
        using RowVector = RMatrixBlock<Derived, 1, Dynamic>;
        using ColVector = RMatrixBlock<Derived, Dynamic, 1>;
    public:
        /* Operations */
        template<class OtherDerived>
        void assignTo(LValueMatrix<OtherDerived>& target) const;
        [[nodiscard]] inline RowVector row(size_t r);
        [[nodiscard]] inline const RowVector row(size_t r) const;
        [[nodiscard]] inline ColVector col(size_t c);
        [[nodiscard]] inline const ColVector col(size_t c) const;
        [[nodiscard]] inline RMatrixBlock<Derived> rows(size_t fromRow, size_t rowCount);
        [[nodiscard]] inline const RMatrixBlock<Derived> rows(size_t fromRow, size_t rowCount) const;
        [[nodiscard]] inline RMatrixBlock<Derived> cols(size_t fromCol, size_t colCount);
        [[nodiscard]] inline const RMatrixBlock<Derived> cols(size_t fromCol, size_t colCount) const;
        [[nodiscard]] inline RMatrixBlock<Derived> topRows(size_t to);
        [[nodiscard]] inline const RMatrixBlock<Derived> topRows(size_t to) const;
        [[nodiscard]] inline RMatrixBlock<Derived> bottomRows(size_t from);
        [[nodiscard]] inline const RMatrixBlock<Derived> bottomRows(size_t from) const;
        [[nodiscard]] inline RMatrixBlock<Derived> leftCols(size_t to);
        [[nodiscard]] inline const RMatrixBlock<Derived> leftCols(size_t to) const;
        [[nodiscard]] inline RMatrixBlock<Derived> rightCols(size_t from);
        [[nodiscard]] inline const RMatrixBlock<Derived> rightCols(size_t from) const;
        [[nodiscard]] inline RMatrixBlock<Derived> topLeftCorner(size_t toRow, size_t toCol);
        [[nodiscard]] inline const RMatrixBlock<Derived> topLeftCorner(size_t toRow, size_t toCol) const;
        [[nodiscard]] inline RMatrixBlock<Derived> topLeftCorner(size_t to);
        [[nodiscard]] inline const RMatrixBlock<Derived> topLeftCorner(size_t to) const;
        [[nodiscard]] inline RMatrixBlock<Derived> topRightCorner(size_t toRow, size_t fromCol);
        [[nodiscard]] inline const RMatrixBlock<Derived> topRightCorner(size_t toRow, size_t fromCol) const;
        [[nodiscard]] inline RMatrixBlock<Derived> bottomLeftCorner(size_t fromRow, size_t toCol);
        [[nodiscard]] inline const RMatrixBlock<Derived> bottomLeftCorner(size_t fromRow, size_t toCol) const;
        [[nodiscard]] inline RMatrixBlock<Derived> bottomRightCorner(size_t fromRow, size_t fromCol);
        [[nodiscard]] inline const RMatrixBlock<Derived> bottomRightCorner(size_t fromRow, size_t fromCol) const;
        [[nodiscard]] inline RMatrixBlock<Derived> bottomRightCorner(size_t from);
        [[nodiscard]] inline const RMatrixBlock<Derived> bottomRightCorner(size_t from) const;
        [[nodiscard]] inline RMatrixBlock<Derived> block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount);
        [[nodiscard]] inline const RMatrixBlock<Derived> block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) const;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return Base::getDerived().calc(row, col); }
        [[nodiscard]] size_t getRow() const noexcept { return Base::getDerived().getRow(); }
        [[nodiscard]] size_t getColumn() const noexcept { return Base::getDerived().getColumn(); }
        [[nodiscard]] ScalarType trace() const;
        [[nodiscard]] Transpose<Derived> transpose() const noexcept;
        [[nodiscard]] Conjugate<Derived> conjugate() const noexcept;
    };

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
        return Transpose<Derived>(this->getDerived());
    }

    template<class Derived>
    Conjugate<Derived> RValueMatrix<Derived>::conjugate() const noexcept {
        return Conjugate<Derived>(this->getDerived());
    }

    template<class Derived>
    std::ostream& operator<<(std::ostream& os, const RValueMatrix<Derived>& m);
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class Derived>
    Derived reciprocal(const RValueMatrix<Derived>& m);
    
    template<class Derived>
    Derived sqrt(const RValueMatrix<Derived>& m);
    
    template<class Derived>
    Derived factorial(const RValueMatrix<Derived>& m);
    
    template<class Derived>
    Derived ln(const RValueMatrix<Derived>& m);
    
    template<class Derived>
    Derived log(const RValueMatrix<Derived>& m, const MultiScalar& a);

    template<class Derived>
    Derived pow(const RValueMatrix<Derived>& m, const MultiScalar& a);

    template<class Derived>
    Derived tan(const RValueMatrix<Derived>& m);
    
    template<class Derived>
    Derived sec(const RValueMatrix<Derived>& m);
    
    template<class Derived>
    Derived csc(const RValueMatrix<Derived>& m);
    
    template<class Derived>
    Derived cot(const RValueMatrix<Derived>& m);
    
    template<class Derived>
    Derived arccos(const RValueMatrix<Derived>& m);
    
    template<class Derived>
    Derived arcsin(const RValueMatrix<Derived>& m);
    
    template<class Derived>
    Derived arctan(const RValueMatrix<Derived>& m);
    
    template<class Derived>
    Derived arcsec(const RValueMatrix<Derived>& m);
    
    template<class Derived>
    Derived arccsc(const RValueMatrix<Derived>& m);
    
    template<class Derived>
    Derived arccot(const RValueMatrix<Derived>& m);
    
    template<class Derived>
    Derived cosh(const RValueMatrix<Derived>& m);
    
    template<class Derived>
    Derived sinh(const RValueMatrix<Derived>& m);
    
    template<class Derived>
    Derived tanh(const RValueMatrix<Derived>& m);
    
    template<class Derived>
    Derived sech(const RValueMatrix<Derived>& m);
    
    template<class Derived>
    Derived csch(const RValueMatrix<Derived>& m);
    
    template<class Derived>
    Derived coth(const RValueMatrix<Derived>& m);
    
    template<class Derived>
    Derived arccosh(const RValueMatrix<Derived>& m);
    
    template<class Derived>
    Derived arcsinh(const RValueMatrix<Derived>& m);
    
    template<class Derived>
    Derived arctanh(const RValueMatrix<Derived>& m);
    
    template<class Derived>
    Derived arcsech(const RValueMatrix<Derived>& m);
    
    template<class Derived>
    Derived arccsch(const RValueMatrix<Derived>& m);
    
    template<class Derived>
    Derived arccoth(const RValueMatrix<Derived>& m);
}

#include "RValueMatrixImpl.h"