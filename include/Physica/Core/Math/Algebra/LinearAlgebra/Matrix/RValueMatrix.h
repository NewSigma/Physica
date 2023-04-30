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

#include "MatrixImpl/RMatrixBlock.h"

namespace Physica::Core {
    template<class Derived> class LValueMatrix;
    template<class MatrixType> class Transpose;
    template<class MatrixType> class Conjugate;
    template<class MatrixType> class Flatten;
    template<class MatrixType, bool isLValueMatrix> class DiagVector;
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
        constexpr static int Option = Internal::Traits<Derived>::Option;
        constexpr static size_t RowAtCompile = Internal::Traits<Derived>::RowAtCompile;
        constexpr static size_t ColumnAtCompile = Internal::Traits<Derived>::ColumnAtCompile;
        constexpr static size_t MaxRowAtCompile = Internal::Traits<Derived>::MaxRowAtCompile;
        constexpr static size_t MaxColumnAtCompile = Internal::Traits<Derived>::MaxColumnAtCompile;
        constexpr static size_t SizeAtCompile = Internal::Traits<Derived>::SizeAtCompile;
        constexpr static size_t MaxSizeAtCompile = Internal::Traits<Derived>::MaxSizeAtCompile;

        constexpr static bool isColumnMatrix = MatrixOption::isColumnMatrix<Derived>();
        constexpr static bool isRowMatrix = MatrixOption::isRowMatrix<Derived>();
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
        [[nodiscard]] inline DiagVector<Derived, false> diag();
        [[nodiscard]] inline const DiagVector<Derived, false> diag() const;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return Base::getDerived().calc(row, col); }
        [[nodiscard]] ScalarType calcFromMajorMinor(size_t row, size_t col) const;
        [[nodiscard]] size_t getRow() const noexcept { return Base::getDerived().getRow(); }
        [[nodiscard]] size_t getColumn() const noexcept { return Base::getDerived().getColumn(); }
        [[nodiscard]] size_t getMaxMajor() const noexcept { return MatrixOption::getMaxMajor<Derived>(Base::getDerived()); }
        [[nodiscard]] size_t getMaxMinor() const noexcept { return MatrixOption::getMaxMinor<Derived>(Base::getDerived()); }
        [[nodiscard]] ScalarType trace() const;
        [[nodiscard]] Transpose<Derived> transpose() const noexcept;
        [[nodiscard]] Conjugate<Derived> conjugate() const noexcept;
        [[nodiscard]] Flatten<Derived> flatten() const noexcept;
        [[nodiscard]] ScalarType sum() const { return Base::getDerived().sum(); }
        /* Static members */
        [[nodiscard]] static size_t rowFromMajorMinor(size_t major, size_t minor) noexcept { return MatrixOption::rowFromMajorMinor<Derived>(major, minor); }
        [[nodiscard]] static size_t columnFromMajorMinor(size_t major, size_t minor) noexcept { return MatrixOption::columnFromMajorMinor<Derived>(major, minor); }
    };

    template<class MatrixType, class MatrixType2>
    bool matrixNear(const RValueMatrix<MatrixType>& m1, const RValueMatrix<MatrixType2>& m2, double precision);

    template<class Derived>
    std::ostream& operator<<(std::ostream& os, const RValueMatrix<Derived>& m);
}

#include "MatrixImpl/RValueMatrixImpl.h"
#ifdef PHYSICA_CUDA
    #include "RValueMatrix.cuh"
#endif
