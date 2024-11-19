/*
 * Copyright 2021-2024 Weibo He.
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
#include "LValueMatrixImpl/LMatrixBlock.h"

namespace Physica::Core {
    template<class MatrixType> class InverseMatrix;
    template<class MatrixType> class LValueFlatten;
    /**
     * \class LValueMatrix is base class of matrixes that can be assigned to \class LValueMatrix
     * and other matrixes can be assigned to this class.
     * In other words, you can take the address of elements in a LValueMatrix.
     */
    template<class Derived>
    class LValueMatrix : public RValueMatrix<Derived> {
        using Base = RValueMatrix<Derived>;
        using This = LValueMatrix<Derived>;
        using RowVector = LMatrixBlock<Derived, 1, Dynamic>;
        using ColVector = LMatrixBlock<Derived, Dynamic, 1>;
    public:
        using typename Base::ScalarType;
        using typename Base::RealType;
        using Base::RowAtCompile;
        using Base::ColAtCompile;
        using Base::isComplex;
        using Base::isForwardDiff;
        using Base::isReverseDiff;
    protected:
        using PtrTy = typename ScalarType::PtrTy;
        using ConstPtrTy = typename ScalarType::ConstPtrTy;
        using RefTy = typename ScalarType::RefTy;
        using ConstRefTy = typename ScalarType::ConstRefTy;
    public:
        ~LValueMatrix() = default;
        /* Operators */
        inline This& operator=(const This& m);
        inline This& operator=(This&& m);

        Derived& operator=(const ScalarType& s);
        void operator+=(const ScalarType& s) { Base::getDerived() = Base::getDerived() + s; }
        void operator-=(const ScalarType& s) { Base::getDerived() = Base::getDerived() - s; }
        void operator*=(const ScalarType& s) { Base::getDerived() = Base::getDerived() * s; }
        void operator/=(const ScalarType& s) { Base::getDerived() = Base::getDerived() / s; }

        template<Matrix M> Derived& operator=(const M& m);
        template<Matrix M> void operator+=(const M& m) { Base::getDerived() = Base::getDerived() + m; }
        template<Matrix M> void operator-=(const M& m) { Base::getDerived() = Base::getDerived() - m; }
        template<Matrix M> void operator*=(const M& m) { Base::getDerived() = Derived(Base::getDerived() * m); }
        [[nodiscard]] inline RefTy operator()(size_t row, size_t col);
        [[nodiscard]] inline ConstRefTy operator()(size_t row, size_t col) const;
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
        void colReduce(size_t c1, size_t c2, size_t elementIndex);
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
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const;
        [[nodiscard]] __host__ __device__ inline PtrTy data_ptr(size_t row, size_t col) noexcept;
        [[nodiscard]] __host__ __device__ inline ConstPtrTy data_ptr(size_t row, size_t col) const noexcept;
        [[nodiscard]] inline RefTy refFromMajorMinor(size_t major, size_t minor);
        [[nodiscard]] inline ConstRefTy refFromMajorMinor(size_t major, size_t minor) const;
        [[nodiscard]] LValueFlatten<Derived> flatten();
        [[nodiscard]] const LValueFlatten<Derived> flatten() const;
        /* Setters */
        void toUnitMatrix();
    protected:
        LValueMatrix() = default;
        LValueMatrix(const This&) = default;
        LValueMatrix(This&&) noexcept = default;
    };
}

namespace Physica {
    template<class Derived>
    class Traits<Core::LValueMatrix<Derived>> : public Traits<Derived> {};
}

#include "LValueMatrixImpl/LValueMatrixImpl.h"
#include "InverseMatrix.h"
