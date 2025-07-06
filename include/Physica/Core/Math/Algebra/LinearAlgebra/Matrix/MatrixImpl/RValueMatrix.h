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

namespace Physica {
    template<class Derived> class LValueMatrix;
    template<class Derived> class ContinuousMatrix;
    template<class, bool ReduceCol> class MatrixSum;
    template<class MatrixType, bool isLValueMatrix> class DiagVector;
    template<class> class Inverse;
    template<class> class Transpose;
    template<class> class Conjugate;
    template<class> class Hermite;
    template<class> class FlattenR;
    template<class> class TrigUpper;
    template<class> class TrigLower;

    template<class> class RealMatrix;
    template<class> class ImagMatrix;
    template<class> class SquaredNormMatrix;
    template<class> class NormMatrix;
    template<class> class ValueMatrix;
    template<class, int GradOrder> class GradMatrix;
    template<Scalar, bool Pivot> class QRDecomp;
    template<Matrix, Matrix> class GEMM;
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
        constexpr static int Option = Traits<Derived>::Option;
        constexpr static size_t RowAtCompile = Traits<Derived>::RowAtCompile;
        constexpr static size_t ColAtCompile = Traits<Derived>::ColAtCompile;
        constexpr static size_t SizeAtCompile = Traits<Derived>::SizeAtCompile;
        constexpr static bool isForwardDiff = ScalarType::isForwardDiff;
        constexpr static bool isReverseDiff = ScalarType::isReverseDiff;
        constexpr static bool isDiffable = ScalarType::isDiffable;
        constexpr static bool isComplex = ScalarType::isComplex;
    protected:
        using T = ScalarType;
        using Tr = T::RealType;
        using Tv = T::ValueType;
        using Trv = Tr::ValueType;
        using Tc = T::ComplexType;
    private:
        using RealsRtnTy = std::conditional<isComplex, RealMatrix<Derived>, Derived&>::type;
        using ValuesRtnTy = std::conditional<isDiffable, ValueMatrix<Derived>, Derived&>::type;
    public:
        ~RValueMatrix() = default;
        /* Operators */
        template<Vector V>
        [[nodiscard]] auto operator*(V&& v) const& noexcept requires(RowAtCompile != 1 && !CUDA<V>);
        template<Vector V>
        [[nodiscard]] auto operator*(V&& v) && noexcept requires(RowAtCompile != 1 && !CUDA<V>);
        template<Vector V>
        [[nodiscard]] auto operator*(const V& v) const noexcept requires(RowAtCompile == 1 && !CUDA<V>);
        template<Matrix M>
        [[nodiscard]] auto operator*(const M& mat) const noexcept requires(((ColAtCompile != 1 && M::ColAtCompile != 1) || (ColAtCompile == 1 && M::ColAtCompile == 1)) && !CUDA<M>);
        /* Operations */
        void assign(Matrix auto& target) const;
        void assign_add(Matrix auto& target) const;

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
        [[nodiscard]] auto triu() noexcept;
        [[nodiscard]] const auto triu() const noexcept;
        [[nodiscard]] auto tril() noexcept;
        [[nodiscard]] const auto tril() const noexcept;

        [[nodiscard]] auto calc(size_t row, size_t col) const { return Base::getDerived().calc(row, col); }
        [[nodiscard]] auto calc_value(size_t row, size_t col) const { return Base::getDerived().calc_value(row, col); }
        [[nodiscard]] inline auto calcFromMajorMinor(size_t major, size_t minor) const;
        void reverse(const Matrix auto& y, const Matrix auto& grad) const noexcept requires(isReverseDiff);

        [[nodiscard]] Tr norm1() const;
        template<ExecutePolicy P = Sequential>
        [[nodiscard]] Tr norm1_power(unsigned int maxIteration) const;
        [[nodiscard]] Tr normF() const;
        [[nodiscard]] Tr normInf() const;
        [[nodiscard]] T cond2() const;

        [[nodiscard]] T max() const;
        [[nodiscard]] T min() const;
        [[nodiscard]] CoDiff<T> sum() const;
        [[nodiscard]] auto sum_cols() const;
        [[nodiscard]] T mean() const;
        [[nodiscard]] T trace() const;
        [[nodiscard]] CoDiff<T> lnSumExp() const;
        [[nodiscard]] CoDiff<T> det() const;
        [[nodiscard]] T lnAbsDet() const;

        [[nodiscard]] auto format() const noexcept;
        [[nodiscard]] auto inverse() const noexcept;
        [[nodiscard]] auto transpose() const noexcept;
        [[nodiscard]] auto conjugate() const noexcept;
        [[nodiscard]] auto hermite() const noexcept;
        [[nodiscard]] auto flatten() const noexcept;

        [[nodiscard]] RealsRtnTy reals() const noexcept;
        [[nodiscard]] auto imags() const noexcept;
        [[nodiscard]] auto squaredNorms() const noexcept;
        [[nodiscard]] auto norms() const noexcept;
        [[nodiscard]] ValuesRtnTy values() const noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] auto grads() const noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return Base::getDerived().getRow(); }
        [[nodiscard]] size_t getCol() const noexcept { return Base::getDerived().getCol(); }
        [[nodiscard]] size_t getMaxMajor() const noexcept { return MatrixOption::getMaxMajor<Derived>(Base::getDerived()); }
        [[nodiscard]] size_t getMaxMinor() const noexcept { return MatrixOption::getMaxMinor<Derived>(Base::getDerived()); }

        [[nodiscard]] bool isSquare() const noexcept;
        [[nodiscard]] bool isSymm() const noexcept;
        [[nodiscard]] bool isFinite() const noexcept;
        /* Static members */
        [[nodiscard]] static size_t rowFromMajorMinor(size_t major, size_t minor) noexcept { return MatrixOption::rowFromMajorMinor<Derived>(major, minor); }
        [[nodiscard]] static size_t colFromMajorMinor(size_t major, size_t minor) noexcept { return MatrixOption::colFromMajorMinor<Derived>(major, minor); }
        template<Matrix M>
        __host__ __device__ static void assign_check(const M& target) noexcept;
        template<Matrix M>
        __host__ __device__ constexpr static bool matmul_check() noexcept;
    protected:
        RValueMatrix() = default;
        RValueMatrix(const This&) = default;
        RValueMatrix(This&&) noexcept = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        /* Operations */
        template<int GradOrder>
        auto grads_impl() const noexcept;
    };

    bool matrixNear(const Matrix auto& m1, const Matrix auto& m2, double precision);

    template<Matrix T>
    bool operator==(const T& m1, const T& m2) {
        if (m1.getRow() != m2.getRow())
            return false;
        if (m1.getCol() != m2.getCol())
            return false;
        for (size_t major = 0; major < m1.getMaxMajor(); ++major)
            for (size_t minor = 0; minor < m1.getMaxMinor(); ++minor)
                if (m1.calcFromMajorMinor(major, minor) != m2.calcFromMajorMinor(major, minor))
                    return false;
        return true;
    }

    template<Matrix T>
    bool operator!=(const T& m1, const T& m2) { return !(m1 == m2); }
}

namespace Physica {
    template<class T>
    class Traits<RValueMatrix<T>> : public Traits<T> {
    public:
        using Derived = T;
    };
}

#include "RValueMatrixImpl/RValueMatrixImpl.h"
#include "MatrixProduct/GEMM.h"
#include "MatrixProduct/GEMV.h"
#include "MatrixProduct/GEVM.h"
#include "RValueMatrixImpl/MatrixNorm.h"
#include "MatrixExpr.h"
#include "DiagVector.h"
