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

#include "RValueMatrix.h"
#include "MatrixBlock.h"

namespace Physica::Core {
    template<class MatrixType> class InverseMatrix;

    namespace Internal {
        template<class Derived>
        class Traits<LValueMatrix<Derived>> : public Traits<Derived> {};
    }
    /**
     * \class LValueMatrix is base class of matrixes that can be assigned to \class LValueMatrix
     * and other matrixes can be assigned to this class.
     * In other words, you can take the address of elements in the matrix.
     */
    template<class Derived>
    class LValueMatrix : public RValueMatrix<Derived> {
    public:
        using Base = RValueMatrix<Derived>;
        using typename Base::ScalarType;
        using RowVector = MatrixBlock<Derived, 1, Dynamic>;
        using ColVector = MatrixBlock<Derived, Dynamic, 1>;
        constexpr static int MatrixOption = Internal::Traits<Derived>::MatrixOption; //It is declared here because MatrixOption makes no sence to a RValueMatrix
        constexpr static bool isColumnMatrix = DenseMatrixOption::isColumnMatrix<Derived>();
        constexpr static bool isRowMatrix = DenseMatrixOption::isRowMatrix<Derived>();
    public:
        ~LValueMatrix() = default;
        /* Operators */
        LValueMatrix& operator=(const LValueMatrix& m) = delete;
        LValueMatrix& operator=(LValueMatrix&& m) = delete;
        template<class OtherMatrix>
        Derived& operator=(const RValueMatrix<OtherMatrix>& m);
        template<ScalarOption option, bool errorTrack>
        Derived& operator=(const Scalar<option, errorTrack>& s); //No need for complex assignment, assign a complex number to real matrix is undefined
        [[nodiscard]] ScalarType& operator()(size_t row, size_t column) { return Base::getDerived()(row, column); }
        [[nodiscard]] const ScalarType& operator()(size_t row, size_t column) const { return Base::getDerived()(row, column); }
        /* Operations */
        [[nodiscard]] inline RowVector row(size_t r);
        [[nodiscard]] inline const RowVector row(size_t r) const;
        [[nodiscard]] inline ColVector col(size_t c);
        [[nodiscard]] inline const ColVector col(size_t c) const;
        [[nodiscard]] inline MatrixBlock<Derived> rows(size_t fromRow, size_t rowCount);
        [[nodiscard]] inline const MatrixBlock<Derived> rows(size_t fromRow, size_t rowCount) const;
        [[nodiscard]] inline MatrixBlock<Derived> cols(size_t fromCol, size_t colCount);
        [[nodiscard]] inline const MatrixBlock<Derived> cols(size_t fromCol, size_t colCount) const;
        [[nodiscard]] inline MatrixBlock<Derived> topRows(size_t to);
        [[nodiscard]] inline const MatrixBlock<Derived> topRows(size_t to) const;
        [[nodiscard]] inline MatrixBlock<Derived> bottomRows(size_t from);
        [[nodiscard]] inline const MatrixBlock<Derived> bottomRows(size_t from) const;
        [[nodiscard]] inline MatrixBlock<Derived> leftCols(size_t from);
        [[nodiscard]] inline const MatrixBlock<Derived> leftCols(size_t from) const;
        [[nodiscard]] inline MatrixBlock<Derived> rightCols(size_t from);
        [[nodiscard]] inline const MatrixBlock<Derived> rightCols(size_t from) const;
        [[nodiscard]] inline MatrixBlock<Derived> topLeftCorner(size_t toRow, size_t toCol);
        [[nodiscard]] inline const MatrixBlock<Derived> topLeftCorner(size_t toRow, size_t toCol) const;
        [[nodiscard]] inline MatrixBlock<Derived> topLeftCorner(size_t to);
        [[nodiscard]] inline const MatrixBlock<Derived> topLeftCorner(size_t to) const;
        [[nodiscard]] inline MatrixBlock<Derived> topRightCorner(size_t toRow, size_t fromCol);
        [[nodiscard]] inline const MatrixBlock<Derived> topRightCorner(size_t toRow, size_t fromCol) const;
        [[nodiscard]] inline MatrixBlock<Derived> bottomLeftCorner(size_t fromRow, size_t toCol);
        [[nodiscard]] inline const MatrixBlock<Derived> bottomLeftCorner(size_t fromRow, size_t toCol) const;
        [[nodiscard]] inline MatrixBlock<Derived> bottomRightCorner(size_t fromRow, size_t fromCol);
        [[nodiscard]] inline const MatrixBlock<Derived> bottomRightCorner(size_t fromRow, size_t fromCol) const;
        [[nodiscard]] inline MatrixBlock<Derived> bottomRightCorner(size_t from);
        [[nodiscard]] inline const MatrixBlock<Derived> bottomRightCorner(size_t from) const;
        [[nodiscard]] inline MatrixBlock<Derived> block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount);
        [[nodiscard]] inline const MatrixBlock<Derived> block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) const;

        [[nodiscard]] InverseMatrix<Derived> inverse() const noexcept;
        ScalarType determinate() const;
        void rowReduce(size_t r1, size_t r2, size_t elementIndex);
        void rowReduce(size_t r1, size_t r2, const ScalarType& factor);
        void columnReduce(size_t c1, size_t c2, size_t elementIndex);
        void columnReduce(size_t c1, size_t c2, const ScalarType& factor);
        inline void majorReduce(size_t v1, size_t v2, size_t elementIndex);
        inline void majorReduce(size_t v1, size_t v2, const ScalarType& factor);
        void rowMulScalar(size_t r, const ScalarType& factor);
        void columnMulScalar(size_t c, const ScalarType& factor);
        inline void majorMulScalar(size_t v, const ScalarType& factor);
        inline void majorSwap(size_t v1, size_t v2);

        template<class OtherDerived>
        void assignTo(LValueMatrix<OtherDerived>& target) const;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return (*this)(row, col); }
        [[nodiscard]] ScalarType max() const;
        [[nodiscard]] ScalarType min() const;
        [[nodiscard]] ScalarType sum() const;
        [[nodiscard]] inline size_t getMaxMajor() const noexcept;
        [[nodiscard]] inline size_t getMaxMinor() const noexcept;
        [[nodiscard]] ScalarType& getElementFromMajorMinor(size_t major, size_t minor);
        [[nodiscard]] const ScalarType& getElementFromMajorMinor(size_t major, size_t minor) const;
        /* Setters */
        void toUnitMatrix();
        /* Static members */
        [[nodiscard]] inline static size_t rowFromMajorMinor(size_t major, size_t minor) noexcept;
        [[nodiscard]] inline static size_t columnFromMajorMinor(size_t major, size_t minor) noexcept;
    protected:
        LValueMatrix() = default;
        LValueMatrix(const LValueMatrix&) = default;
        LValueMatrix(LValueMatrix&&) noexcept = default;
    };

    template<class Derived, class OtherDerived>
    inline void operator+=(LValueMatrix<Derived>& m1, const RValueMatrix<OtherDerived>& m2) { m1 = m1 + m2; }
    template<class MatrixType, ScalarOption option, bool errorTrack>
    inline void operator+=(LValueMatrix<MatrixType>& m, const Scalar<option, errorTrack>& s) { m = m + s; }
    template<class Derived, class OtherDerived>
    inline void operator-=(LValueMatrix<Derived>& m1, const RValueMatrix<OtherDerived>& m2) { m1 = m1 - m2; }
    template<class MatrixType, ScalarOption option, bool errorTrack>
    inline void operator-=(LValueMatrix<MatrixType>& m, const Scalar<option, errorTrack>& s) { m = m - s; }
    template<class Derived, class OtherDerived>
    inline void operator*=(LValueMatrix<Derived>& m1, const RValueMatrix<OtherDerived>& m2) { m1 = m1 * m2; }
    template<class MatrixType, ScalarOption option, bool errorTrack>
    inline void operator*=(LValueMatrix<MatrixType>& m, const Scalar<option, errorTrack>& s) { m = m * s; }
}

#include "LValueMatrixImpl.h"
