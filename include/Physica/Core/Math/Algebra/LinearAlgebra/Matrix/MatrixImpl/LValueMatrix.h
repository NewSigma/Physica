/*
 * Copyright 2021-2025 Weibo He.
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

namespace Physica {
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
        using BlockType = LMatrixBlock<Derived>;
    public:
        using typename Base::ScalarType;
        using Base::RowAtCompile;
        using Base::ColAtCompile;
        using Base::isComplex;
        using Base::isForwardDiff;
        using Base::isReverseDiff;
    protected:
        using typename Base::Tr;
        using typename Base::Tv;
        using PtrTy = ScalarType::PtrTy;
        using ConstPtrTy = ScalarType::ConstPtrTy;
        using RefTy = ScalarType::RefTy;
        using ConstRefTy = ScalarType::ConstRefTy;
    public:
        ~LValueMatrix() = default;
        /* Operators */
        inline This& operator=(const This& m);
        inline This& operator=(This&& m);

        template<Scalar T> Derived& operator=(const T& x) requires(!isReverseDiff || !ReverseDiff<T>);
        template<Scalar T> void operator+=(const T& x) { Base::getDerived() = Base::getDerived() + x; }
        template<Scalar T> void operator-=(const T& x) { Base::getDerived() = Base::getDerived() - x; }
        template<Scalar T> void operator*=(const T& x) { Base::getDerived() = Base::getDerived() * x; }
        template<Scalar T> void operator/=(const T& x) { Base::getDerived() = Base::getDerived() / x; }

        template<Matrix M> Derived& operator=(const M& m);
        template<Matrix M> void operator+=(const M& m) { Base::getDerived() = Base::getDerived() + m; }
        template<Matrix M> void operator-=(const M& m) { Base::getDerived() = Base::getDerived() - m; }
        template<Matrix M> void operator*=(const M& m) { Base::getDerived() = Derived(Base::getDerived() * m); }
        [[nodiscard]] inline RefTy operator()(size_t row, size_t col);
        [[nodiscard]] inline ConstRefTy operator()(size_t row, size_t col) const;
        /* Operations */
        template<Matrix M>
        void reverse(const M& grad) const noexcept requires(isReverseDiff);

        [[nodiscard]] inline auto row(size_t r) noexcept;
        [[nodiscard]] inline const auto row(size_t r) const noexcept;
        [[nodiscard]] inline auto col(size_t c) noexcept;
        [[nodiscard]] inline const auto col(size_t c) const noexcept;
        [[nodiscard]] inline auto rows(size_t fromRow, size_t rowCount) noexcept;
        [[nodiscard]] inline const auto rows(size_t fromRow, size_t rowCount) const noexcept;
        [[nodiscard]] inline auto topRows(size_t to) noexcept;
        [[nodiscard]] inline const auto topRows(size_t to) const noexcept;
        [[nodiscard]] inline auto bottomRows(size_t from) noexcept;
        [[nodiscard]] inline const auto bottomRows(size_t from) const noexcept;
        [[nodiscard]] inline auto cols(size_t fromCol, size_t colCount) noexcept;
        [[nodiscard]] inline const auto cols(size_t fromCol, size_t colCount) const noexcept;
        [[nodiscard]] inline auto leftCols(size_t to) noexcept;
        [[nodiscard]] inline const auto leftCols(size_t to) const noexcept;
        [[nodiscard]] inline auto rightCols(size_t from) noexcept;
        [[nodiscard]] inline const auto rightCols(size_t from) const noexcept;
        [[nodiscard]] inline auto topLeftCorner(size_t toRow, size_t toCol) noexcept;
        [[nodiscard]] inline const auto topLeftCorner(size_t toRow, size_t toCol) const noexcept;
        [[nodiscard]] inline auto topLeftCorner(size_t to) noexcept;
        [[nodiscard]] inline const auto topLeftCorner(size_t to) const noexcept;
        [[nodiscard]] inline auto topRightCorner(size_t toRow, size_t fromCol) noexcept;
        [[nodiscard]] inline const auto topRightCorner(size_t toRow, size_t fromCol) const noexcept;
        [[nodiscard]] inline auto bottomLeftCorner(size_t fromRow, size_t toCol) noexcept;
        [[nodiscard]] inline const auto bottomLeftCorner(size_t fromRow, size_t toCol) const noexcept;
        [[nodiscard]] inline auto bottomRightCorner(size_t fromRow, size_t fromCol) noexcept;
        [[nodiscard]] inline const auto bottomRightCorner(size_t fromRow, size_t fromCol) const noexcept;
        [[nodiscard]] inline auto bottomRightCorner(size_t from) noexcept;
        [[nodiscard]] inline const auto bottomRightCorner(size_t from) const noexcept;
        [[nodiscard]] inline auto block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept;
        [[nodiscard]] inline const auto block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) const noexcept;
        [[nodiscard]] inline auto diag() noexcept;
        [[nodiscard]] inline const auto diag() const noexcept;

        [[nodiscard]] ConstRefTy calc(size_t row, size_t col) const { return operator()(row, col); }
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const { return calc(row, col).value(); }

        [[nodiscard]] CoDiff<ScalarType> sum() const;

        [[nodiscard]] auto inverse() const noexcept;
        [[nodiscard]] CoDiff<ScalarType> determinate() const;
        void rowReduce(size_t r1, size_t r2, size_t elementIndex);
        void colReduce(size_t c1, size_t c2, size_t elementIndex);
        inline void majorReduce(size_t v1, size_t v2, size_t elementIndex);
        inline void majorReduce(size_t v1, size_t v2, const ScalarType& factor);
        inline void majorMulScalar(size_t v, const ScalarType& factor);
        inline void majorSwap(size_t v1, size_t v2);

        template<RNG R>
        void random_uniform();
        template<RNG R>
        void random_normal();
        template<RNG R, class Distribution>
        void random_any(Distribution& dist);

        template<int GradOrder = 1>
        auto grads() const noexcept;
        /* Getters */
        [[nodiscard]] inline PtrTy data_ptr(size_t row, size_t col) noexcept;
        [[nodiscard]] inline ConstPtrTy data_ptr(size_t row, size_t col) const noexcept;
        [[nodiscard]] inline RefTy refFromMajorMinor(size_t major, size_t minor);
        [[nodiscard]] inline ConstRefTy refFromMajorMinor(size_t major, size_t minor) const;
        [[nodiscard]] auto flatten();
        [[nodiscard]] const auto flatten() const;
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
    class Traits<LValueMatrix<Derived>> : public Traits<Derived> {};
}

#include "LValueMatrixImpl/LValueMatrixImpl.h"
#include "LValueMatrixImpl/ReshapedVector.h"
#include "LValueMatrixImpl/InverseMatrix.h"
