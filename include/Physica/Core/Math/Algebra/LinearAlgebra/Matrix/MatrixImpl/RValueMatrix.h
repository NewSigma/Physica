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

#include "RValueMatrixImpl/RMatrixBlock.h"

namespace Physica::Core {
    template<class Derived> class RValueMatrix;
    template<class Derived> class LValueMatrix;
    template<class Derived> class ContinuousMatrix;
    template<class MatrixType> class Transpose;
    template<class MatrixType> class Conjugate;
    template<class MatrixType> class Hermite;
    template<class MatrixType> class RValueFlatten;
    template<class MatrixType> class FormatedMatrix;
    template<class MatrixType, bool isLValueMatrix> class DiagVector;

    template<class T>
    struct is_matrix {
        constexpr static bool value = std::is_base_of<RValueMatrix<T>, T>::value;
    };
    /**
     * \class RValueMatrix is base class of matrixes that can be assigned to \class LValueMatrix
     * and other matrixes cannot be assigned to this class.
     * In other words, you cannot take the address of elements in a RValueMatrix but can calculate its value.
     */
    template<class Derived>
    class RValueMatrix : public CRTPBase<Derived> {
        static_assert(!std::is_const<Derived>::value, "[Error]: A common mistake, const is unnecessary");
        static_assert(!std::is_volatile<Derived>::value, "[Error]: A common mistake, volatile is unnecessary");
        using Base = CRTPBase<Derived>;
        using RowVector = RMatrixBlock<Derived, 1, Dynamic>;
        using ColVector = RMatrixBlock<Derived, Dynamic, 1>;
        using BlockType = RMatrixBlock<Derived>;
    public:
        using ScalarType = typename Traits<Derived>::ScalarType;
        using PlainScalar = typename ScalarType::PlainScalar;
        using RealType = typename ScalarType::RealType;
        constexpr static int Option = Traits<Derived>::Option;
        constexpr static size_t RowAtCompile = Traits<Derived>::RowAtCompile;
        constexpr static size_t ColumnAtCompile = Traits<Derived>::ColumnAtCompile;
        constexpr static size_t SizeAtCompile = Traits<Derived>::SizeAtCompile;
        constexpr static bool isReverseDiff = ScalarType::isReverseDiff;
        constexpr static bool isComplex = ScalarType::isComplex;
        constexpr static bool isColumnMatrix = MatrixOption::isColumnMatrix<Derived>();
        constexpr static bool isRowMatrix = MatrixOption::isRowMatrix<Derived>();
    public:
        /* Operations */
        template<class OtherDerived>
        void assignTo(LValueMatrix<OtherDerived>& target) const;
        [[nodiscard]] inline RowVector row(size_t r);
        [[nodiscard]] inline const RowVector row(size_t r) const;
        [[nodiscard]] inline ColVector col(size_t c);
        [[nodiscard]] inline const ColVector col(size_t c) const;
        [[nodiscard]] inline BlockType rows(size_t fromRow, size_t rowCount);
        [[nodiscard]] inline const BlockType rows(size_t fromRow, size_t rowCount) const;
        [[nodiscard]] inline BlockType topRows(size_t to);
        [[nodiscard]] inline const BlockType topRows(size_t to) const;
        [[nodiscard]] inline BlockType bottomRows(size_t from);
        [[nodiscard]] inline const BlockType bottomRows(size_t from) const;
        [[nodiscard]] inline BlockType cols(size_t fromCol, size_t colCount);
        [[nodiscard]] inline const BlockType cols(size_t fromCol, size_t colCount) const;
        [[nodiscard]] inline BlockType leftCols(size_t to);
        [[nodiscard]] inline const BlockType leftCols(size_t to) const;
        [[nodiscard]] inline BlockType rightCols(size_t from);
        [[nodiscard]] inline const BlockType rightCols(size_t from) const;
        [[nodiscard]] inline BlockType topLeftCorner(size_t toRow, size_t toCol);
        [[nodiscard]] inline const BlockType topLeftCorner(size_t toRow, size_t toCol) const;
        [[nodiscard]] inline BlockType topLeftCorner(size_t to);
        [[nodiscard]] inline const BlockType topLeftCorner(size_t to) const;
        [[nodiscard]] inline BlockType topRightCorner(size_t toRow, size_t fromCol);
        [[nodiscard]] inline const BlockType topRightCorner(size_t toRow, size_t fromCol) const;
        [[nodiscard]] inline BlockType bottomLeftCorner(size_t fromRow, size_t toCol);
        [[nodiscard]] inline const BlockType bottomLeftCorner(size_t fromRow, size_t toCol) const;
        [[nodiscard]] inline BlockType bottomRightCorner(size_t fromRow, size_t fromCol);
        [[nodiscard]] inline const BlockType bottomRightCorner(size_t fromRow, size_t fromCol) const;
        [[nodiscard]] inline BlockType bottomRightCorner(size_t from);
        [[nodiscard]] inline const BlockType bottomRightCorner(size_t from) const;
        [[nodiscard]] inline BlockType block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount);
        [[nodiscard]] inline const BlockType block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) const;
        [[nodiscard]] inline DiagVector<Derived, false> diag();
        [[nodiscard]] inline const DiagVector<Derived, false> diag() const;

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return Base::getDerived().calc(row, col); }
        [[nodiscard]] inline ScalarType calcFromMajorMinor(size_t major, size_t minor) const;
        [[nodiscard]] inline FormatedMatrix<Derived> format() const;

        [[nodiscard]] RealType norm1() const;
        template<class Executor = SequentialExecutor>
        [[nodiscard]] RealType norm1_power(unsigned int maxIteration) const;

        [[nodiscard]] ScalarType max() const;
        [[nodiscard]] ScalarType min() const;
        [[nodiscard]] ScalarType trace() const;
        [[nodiscard]] inline Transpose<Derived> transpose() const noexcept;
        [[nodiscard]] inline Conjugate<Derived> conjugate() const noexcept;
        [[nodiscard]] inline Hermite<Derived> hermite() const noexcept;
        [[nodiscard]] RValueFlatten<Derived> flatten() const noexcept;
        [[nodiscard]] ScalarType sum() const;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return Base::getDerived().getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const noexcept { return Base::getDerived().getColumn(); }
        [[nodiscard]] __host__ __device__ size_t getMaxMajor() const noexcept { return MatrixOption::getMaxMajor<Derived>(Base::getDerived()); }
        [[nodiscard]] __host__ __device__ size_t getMaxMinor() const noexcept { return MatrixOption::getMaxMinor<Derived>(Base::getDerived()); }
        /* Static members */
        [[nodiscard]] __host__ __device__ static size_t rowFromMajorMinor(size_t major, size_t minor) noexcept { return MatrixOption::rowFromMajorMinor<Derived>(major, minor); }
        [[nodiscard]] __host__ __device__ static size_t columnFromMajorMinor(size_t major, size_t minor) noexcept { return MatrixOption::columnFromMajorMinor<Derived>(major, minor); }
    };

    template<class MatrixType, class MatrixType2>
    bool matrixNear(const RValueMatrix<MatrixType>& m1, const RValueMatrix<MatrixType2>& m2, double precision);

    template<class Derived>
    bool operator==(const RValueMatrix<Derived>& m1, const RValueMatrix<Derived>& m2);

    template<class Derived>
    inline bool operator!=(const RValueMatrix<Derived>& m1, const RValueMatrix<Derived>& m2) { return !(m1 == m2); }
}

#include "RValueMatrixImpl/RValueMatrixImpl.h"
#include "RValueMatrixImpl/RValueFlatten.h"
#include "Transpose.h"
#include "Conjugate.h"
#include "Hermite.h"
#include "DiagVector.h"
#include "MatrixNorm.h"
#include "MatrixExpr.h"
#include "MatrixConvert.h"
#include "MatrixProduct/GEMM.h"
#include "MatrixProduct/GEMV.h"
#include "MatrixProduct/GEVM.h"
