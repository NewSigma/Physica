/*
 * Copyright 2021-2026 Weibo He.
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
    template<Matrix> class DiagVectorL;
    template<Matrix M> class MinorDiagL;
    template<class> class FlattenL;
    /**
     * \class LValueMatrix is base class of matrixes that can be assigned to \class LValueMatrix
     * and other matrixes can be assigned to this class.
     * In other words, you can take the address of elements in a LValueMatrix.
     */
    template<class Derived>
    class LValueMatrix : public RValueMatrix<Derived> {
        using Base = RValueMatrix<Derived>;
        using This = LValueMatrix<Derived>;
    public:
        using typename Base::ScalarType;
        using Base::RowAtCompile;
        using Base::ColAtCompile;
        using Base::isComplex;
        using Base::isForwardDiff;
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tr;
        using typename Base::Tv;
        using typename Base::Trv;
    public:
        ~LValueMatrix() = default;
        /* Operators */
        This& operator=(const This& m) = delete;
        This& operator=(This&& m) noexcept = delete;

        auto operator=(Scalar auto x) noexcept -> Derived&;
        void operator+=(Scalar auto x) noexcept;
        void operator-=(Scalar auto x) noexcept;
        void operator*=(Scalar auto x) noexcept;
        void operator/=(Scalar auto x) noexcept;

        auto operator=(const Matrix auto& m) -> Derived&;
        void operator+=(const Matrix auto& m);
        void operator-=(const Matrix auto& m);
        void operator*=(const Matrix auto& m) { Base::getDerived() = Derived(Base::getDerived() * m); }

        [[nodiscard]] decltype(auto) operator[](this auto&&, size_t row, size_t col);
        /* Operations */
        [[nodiscard]] decltype(auto) calc(size_t row, size_t col) const { return operator[](row, col); }
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const { return calc(row, col).value(); }

        [[nodiscard]] size_t pivotPartial(size_t column) const noexcept;
        [[nodiscard]] size_t pivotComplete(size_t column) const noexcept;

        [[nodiscard]] CoDiff<ScalarType> sum() const;

        void toNextMean(size_t lastNumSample, const Matrix auto& sample) noexcept;
        void toNextVariance(Derived& mean, size_t lastNumSample, const Matrix auto& sample) noexcept;
        void reverse(const auto& grad) const noexcept;

        [[nodiscard]] auto row(this auto&&, size_t r) noexcept;
        [[nodiscard]] auto col(this auto&&, size_t c) noexcept;
        [[nodiscard]] auto rows(this auto&&, size_t fromRow, size_t rowCount) noexcept;
        [[nodiscard]] auto topRows(this auto&&, size_t to) noexcept;
        [[nodiscard]] auto bottomRows(this auto&&, size_t from) noexcept;
        [[nodiscard]] auto cols(this auto&&, size_t fromCol, size_t colCount) noexcept;
        [[nodiscard]] auto leftCols(this auto&&, size_t to) noexcept;
        [[nodiscard]] auto rightCols(this auto&&, size_t from) noexcept;
        [[nodiscard]] auto topLeftCorner(this auto&&, size_t toRow, size_t toCol) noexcept;
        [[nodiscard]] auto topLeftCorner(this auto&&, size_t to) noexcept;
        [[nodiscard]] auto topRightCorner(this auto&&, size_t toRow, size_t fromCol) noexcept;
        [[nodiscard]] auto bottomLeftCorner(this auto&&, size_t fromRow, size_t toCol) noexcept;
        [[nodiscard]] auto bottomRightCorner(this auto&&, size_t fromRow, size_t fromCol) noexcept;
        [[nodiscard]] auto bottomRightCorner(this auto&&, size_t from) noexcept;
        [[nodiscard]] auto block(this auto&&, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept;
        [[nodiscard]] auto diag(this auto&&) noexcept;
        [[nodiscard]] auto diag(this auto&&, ssize_t shift) noexcept;

        void swap_row(size_t r1, size_t r2);
        void swap_col(size_t c1, size_t c2);
        void rowReduce(size_t r1, size_t r2, size_t elementIndex);
        void colReduce(size_t c1, size_t c2, size_t elementIndex);
        void majorReduce(size_t v1, size_t v2, size_t elementIndex);
        void majorReduce(size_t v1, size_t v2, const ScalarType& factor);
        void majorMulScalar(size_t v, const ScalarType& factor);
        void majorSwap(size_t v1, size_t v2);

        [[nodiscard]] auto flatten(this auto&&);

        void zeros() noexcept;
        void clamp_min(Tv minimum);
        void clamp_max(Tv maximum);
        template<RNG R>
        void random_uniform();
        template<RNG R>
        void random_normal();
        template<RNG R>
        void random_any(auto& distribution);
        void toIdentity();
        [[nodiscard]] VectorND<T> balance();

        template<int GradOrder = 1>
        auto grads() const noexcept;
        /* Getters */
        [[nodiscard]] auto data_ptr(this auto&&, size_t row, size_t col) noexcept;
        [[nodiscard]] decltype(auto) refFromMajorMinor(this auto&&, size_t major, size_t minor) noexcept;
    protected:
        LValueMatrix() = default;
        LValueMatrix(const This&) = default;
        LValueMatrix(This&&) noexcept = default;
        /* Operations */
        void assert_balance() const noexcept;
    };
}

namespace Physica {
    template<class Derived>
    class Traits<LValueMatrix<Derived>> : public Traits<Derived> {};
}

#include "LValueMatrixImpl/LValueMatrixImpl.h"
#include "LValueMatrixImpl/MinorDiag.h"
#include "LValueMatrixImpl/ReshapedVector.h"
#include "LValueMatrixImpl/DiagVector.h"
