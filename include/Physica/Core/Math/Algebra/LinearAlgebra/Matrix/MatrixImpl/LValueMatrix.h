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

#include "RValueMatrix.h"
#include "LMatrixBlock.h"

namespace Physica::Core {
    template<class MatrixType> class InverseMatrix;
    template<class MatrixType> class LValueFlatten;

    namespace Internal {
        template<class Derived>
        class Traits<LValueMatrix<Derived>> : public Traits<Derived> {};
    }
    /**
     * \class LValueMatrix is base class of matrixes that can be assigned to \class LValueMatrix
     * and other matrixes can be assigned to this class.
     * In other words, you can take the address of elements in a LValueMatrix.
     */
    template<class Derived>
    class LValueMatrix : public RValueMatrix<Derived> {
    public:
        using Base = RValueMatrix<Derived>;
        using typename Base::ScalarType;
        using Base::RowAtCompile;
        using Base::ColumnAtCompile;
        using Base::isComplex;
        using Base::isColumnMatrix;
        using Base::isRowMatrix;
        using RowVector = LMatrixBlock<Derived, 1, Dynamic>;
        using ColVector = LMatrixBlock<Derived, Dynamic, 1>;
    public:
        ~LValueMatrix() = default;
        /* Operators */
        LValueMatrix& operator=(const LValueMatrix& m) = delete;
        LValueMatrix& operator=(LValueMatrix&& m) = delete;
        template<class OtherMatrix> Derived& operator=(const RValueMatrix<OtherMatrix>& m);
        Derived& operator=(const ScalarType& s);
        [[nodiscard]] ScalarType& operator()(size_t row, size_t column) { return *data_ptr(row, column); }
        [[nodiscard]] const ScalarType& operator()(size_t row, size_t column) const { return *data_ptr(row, column); }
        void operator+=(const ScalarType& s) { (*this) = (*this) + s; }
        void operator-=(const ScalarType& s) { (*this) = (*this) - s; }
        void operator*=(const ScalarType& s) { (*this) = (*this) * s; }
        void operator/=(const ScalarType& s) { (*this) = (*this) / s; }
        template<class OtherDerived> void operator+=(const RValueMatrix<OtherDerived>& m) { (*this) = (*this) + m; }
        template<class OtherDerived> void operator-=(const RValueMatrix<OtherDerived>& m) { (*this) = (*this) - m; }
        template<class OtherDerived> void operator*=(const RValueMatrix<OtherDerived>& m) { Base::getDerived() = Derived((*this) * m); }
        /* Operations */
        [[nodiscard]] inline RowVector row(size_t r);
        [[nodiscard]] inline const RowVector row(size_t r) const;
        [[nodiscard]] inline ColVector col(size_t c);
        [[nodiscard]] inline const ColVector col(size_t c) const;
        [[nodiscard]] inline LMatrixBlock<Derived> rows(size_t fromRow, size_t rowCount);
        [[nodiscard]] inline const LMatrixBlock<Derived> rows(size_t fromRow, size_t rowCount) const;
        [[nodiscard]] inline LMatrixBlock<Derived> topRows(size_t to);
        [[nodiscard]] inline const LMatrixBlock<Derived> topRows(size_t to) const;
        [[nodiscard]] inline LMatrixBlock<Derived> bottomRows(size_t from);
        [[nodiscard]] inline const LMatrixBlock<Derived> bottomRows(size_t from) const;
        [[nodiscard]] inline LMatrixBlock<Derived> cols(size_t fromCol, size_t colCount);
        [[nodiscard]] inline const LMatrixBlock<Derived> cols(size_t fromCol, size_t colCount) const;
        [[nodiscard]] inline LMatrixBlock<Derived> leftCols(size_t to);
        [[nodiscard]] inline const LMatrixBlock<Derived> leftCols(size_t to) const;
        [[nodiscard]] inline LMatrixBlock<Derived> rightCols(size_t from);
        [[nodiscard]] inline const LMatrixBlock<Derived> rightCols(size_t from) const;
        [[nodiscard]] inline LMatrixBlock<Derived> topLeftCorner(size_t toRow, size_t toCol);
        [[nodiscard]] inline const LMatrixBlock<Derived> topLeftCorner(size_t toRow, size_t toCol) const;
        [[nodiscard]] inline LMatrixBlock<Derived> topLeftCorner(size_t to);
        [[nodiscard]] inline const LMatrixBlock<Derived> topLeftCorner(size_t to) const;
        [[nodiscard]] inline LMatrixBlock<Derived> topRightCorner(size_t toRow, size_t fromCol);
        [[nodiscard]] inline const LMatrixBlock<Derived> topRightCorner(size_t toRow, size_t fromCol) const;
        [[nodiscard]] inline LMatrixBlock<Derived> bottomLeftCorner(size_t fromRow, size_t toCol);
        [[nodiscard]] inline const LMatrixBlock<Derived> bottomLeftCorner(size_t fromRow, size_t toCol) const;
        [[nodiscard]] inline LMatrixBlock<Derived> bottomRightCorner(size_t fromRow, size_t fromCol);
        [[nodiscard]] inline const LMatrixBlock<Derived> bottomRightCorner(size_t fromRow, size_t fromCol) const;
        [[nodiscard]] inline LMatrixBlock<Derived> bottomRightCorner(size_t from);
        [[nodiscard]] inline const LMatrixBlock<Derived> bottomRightCorner(size_t from) const;
        [[nodiscard]] inline LMatrixBlock<Derived> block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount);
        [[nodiscard]] inline const LMatrixBlock<Derived> block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) const;
        [[nodiscard]] inline DiagVector<Derived, true> diag();
        [[nodiscard]] inline const DiagVector<Derived, true> diag() const;

        [[nodiscard]] InverseMatrix<Derived> inverse() const noexcept;
        ScalarType determinate() const;
        void rowReduce(size_t r1, size_t r2, size_t elementIndex);
        void columnReduce(size_t c1, size_t c2, size_t elementIndex);
        inline void majorReduce(size_t v1, size_t v2, size_t elementIndex);
        inline void majorReduce(size_t v1, size_t v2, const ScalarType& factor);
        inline void majorMulScalar(size_t v, const ScalarType& factor);
        inline void majorSwap(size_t v1, size_t v2);

        template<class RandomGenerator>
        void random_uniform(RandomGenerator& gen);
        template<class RandomGenerator>
        void random_normal(RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        void random_any(Distribution& dist, RandomGenerator& gen);
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return *data_ptr(row, col); }
        [[nodiscard]] __host__ __device__ ScalarType* data_ptr(size_t row, size_t column) { return Base::getDerived().data_ptr(row, column); }
        [[nodiscard]] __host__ __device__ const ScalarType* data_ptr(size_t row, size_t column) const { return Base::getDerived().data_ptr(row, column); }
        [[nodiscard]] inline ScalarType& refFromMajorMinor(size_t major, size_t minor);
        [[nodiscard]] inline const ScalarType& refFromMajorMinor(size_t major, size_t minor) const;
        [[nodiscard]] LValueFlatten<Derived> flatten();
        [[nodiscard]] const LValueFlatten<Derived> flatten() const;
        /* Setters */
        void toUnitMatrix();
    protected:
        LValueMatrix() = default;
        LValueMatrix(const LValueMatrix&) = default;
        LValueMatrix(LValueMatrix&&) noexcept = default;
    };
}

#include "LValueMatrixImpl.h"
