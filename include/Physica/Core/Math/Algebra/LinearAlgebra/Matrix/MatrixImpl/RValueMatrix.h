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

#include "RValueMatrixImpl/RMatrixBlock.h"

namespace Physica::Core {
    template<class Derived> class LValueMatrix;
    template<class Derived> class ContinuousMatrix;
    template<class MatrixType> class Transpose;
    template<class MatrixType> class Conjugate;
    template<class MatrixType> class Hermite;
    template<class MatrixType> class RValueFlatten;
    template<class MatrixType, bool isLValueMatrix> class DiagVector;
    template<class T> class RealMatrix;
    template<class T> class ImagMatrix;
    template<class T> class SquaredNormMatrix;
    template<class T> class NormMatrix;
    template<class T> class ValueMatrix;
    template<class T, int GradOrder> class GradMatrix;
    /**
     * \class RValueMatrix: The base class of all matrixes
     */
    template<class Derived>
    class RValueMatrix : public CRTPBase<RValueMatrix<Derived>> {
        static_assert(!std::is_const<Derived>::value, "[Error]: A common mistake, const is unnecessary");
        static_assert(!std::is_volatile<Derived>::value, "[Error]: A common mistake, volatile is unnecessary");
        using This = RValueMatrix<Derived>;
        using Base = CRTPBase<This>;
        using RowVector = RMatrixBlock<Derived, 1, Dynamic>;
        using ColVector = RMatrixBlock<Derived, Dynamic, 1>;
        using BlockType = RMatrixBlock<Derived>;
    public:
        using ScalarType = Traits<Derived>::ScalarType;
        using ValueType = ScalarType::ValueType;
        using RealType = ScalarType::RealType;
        constexpr static int Option = Traits<Derived>::Option;
        constexpr static size_t RowAtCompile = Traits<Derived>::RowAtCompile;
        constexpr static size_t ColAtCompile = Traits<Derived>::ColAtCompile;
        constexpr static size_t SizeAtCompile = Traits<Derived>::SizeAtCompile;
        constexpr static bool isForwardDiff = ScalarType::isForwardDiff;
        constexpr static bool isReverseDiff = ScalarType::isReverseDiff;
        constexpr static bool isComplex = ScalarType::isComplex;
    public:
        ~RValueMatrix() = default;
        /* Operations */
        template<Matrix M>
        void assign(LValueMatrix<M>& target) const;
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

        [[nodiscard]] auto calc(size_t row, size_t col) const { return Base::getDerived().calc(row, col); }
        [[nodiscard]] auto calc_value(size_t row, size_t col) const { return Base::getDerived().calc_value(row, col); }
        [[nodiscard]] inline auto calcFromMajorMinor(size_t major, size_t minor) const;
        [[nodiscard]] inline auto format() const noexcept;

        [[nodiscard]] RealType norm1() const;
        template<class Executor = SequentialExecutor>
        [[nodiscard]] RealType norm1_power(unsigned int maxIteration) const;
        [[nodiscard]] RealType normInf() const;

        [[nodiscard]] ScalarType max() const;
        [[nodiscard]] ScalarType min() const;
        [[nodiscard]] CoDiff<ScalarType> sum() const;
        [[nodiscard]] ScalarType trace() const;
        [[nodiscard]] inline auto transpose() const noexcept;
        [[nodiscard]] inline auto conjugate() const noexcept;
        [[nodiscard]] inline auto hermite() const noexcept;
        [[nodiscard]] inline auto flatten() const noexcept;

        auto reals() const noexcept;
        auto imags() const noexcept;
        auto squaredNorms() const noexcept;
        auto norms() const noexcept;
        auto values() const noexcept;
        template<int GradOrder = 1>
        auto grads() const noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return Base::getDerived().getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return Base::getDerived().getCol(); }
        [[nodiscard]] __host__ __device__ size_t getMaxMajor() const noexcept { return MatrixOption::getMaxMajor<Derived>(Base::getDerived()); }
        [[nodiscard]] __host__ __device__ size_t getMaxMinor() const noexcept { return MatrixOption::getMaxMinor<Derived>(Base::getDerived()); }

        [[nodiscard]] bool isSymm() const noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ static size_t rowFromMajorMinor(size_t major, size_t minor) noexcept { return MatrixOption::rowFromMajorMinor<Derived>(major, minor); }
        [[nodiscard]] __host__ __device__ static size_t colFromMajorMinor(size_t major, size_t minor) noexcept { return MatrixOption::colFromMajorMinor<Derived>(major, minor); }
    protected:
        RValueMatrix() = default;
        RValueMatrix(const This&) = default;
        RValueMatrix(This&&) noexcept = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
    };

    template<Matrix T1, Matrix T2>
    bool matrixNear(const T1& m1, const T2& m2, double precision);

    template<Matrix T>
    bool operator==(const T& m1, const T& m2);

    template<Matrix T>
    inline bool operator!=(const T& m1, const T& m2) { return !(m1 == m2); }
}

namespace Physica {
    template<class T>
    class Traits<RValueMatrix<T>> : public Traits<T> {
    public:
        using Derived = T;
    };
}

#include "RValueMatrixImpl/RValueMatrixImpl.h"
#include "RValueMatrixImpl/RValueFlatten.h"
#include "RValueMatrixImpl/MatrixConvert.h"
#include "RValueMatrixImpl/ReshapedVector.h"
#include "MatrixProduct/GEMM.h"
#include "MatrixProduct/GEMV.h"
#include "MatrixProduct/GEVM.h"
#include "Transpose.h"
#include "Conjugate.h"
#include "Hermite.h"
#include "DiagVector.h"
#include "MatrixNorm.h"
#include "MatrixExpr.h"
