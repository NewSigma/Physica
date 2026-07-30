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

#include "RValueMatrixImpl/RMatrixBlock.h"

namespace Physica {
    template<class Derived> class LValueMatrix;
    template<class Derived> class CompactMatrix;
    template<class, bool ReduceCol> class MatrixSum;
    template<Matrix> class MainDiag;
    template<Matrix> class OffsetDiag;
    template<class> class Inverse;
    template<Matrix> class PseudoInverse;
    template<class> class Transpose;
    template<class> class Conjugate;
    template<class> class Hermite;
    template<class> class Flatten;
    template<class, bool Upper, bool Unit> class MatrixTrig;

    template<class> class RealMatrix;
    template<class> class ImagMatrix;
    template<class> class SquaredNormMatrix;
    template<class> class NormMatrix;
    template<class> class ValueMatrix;
    template<class, int GradOrder> class GradMatrix;
    template<Scalar, bool Pivot> class DenseLU;
    template<Matrix, Matrix> class GEMM;

    namespace Internal {
        template<Matrix M1, Matrix M2>
        class EnableLAPACK<M1, M2> {
            using U1 = std::remove_cvref<M1>::type;
            using U2 = std::remove_cvref<M2>::type;
            using T1 = U1::ScalarType;
            using T2 = U2::ScalarType;
        public:
            constexpr static bool value = std::same_as<T1, T2>
                                       && !Diffable<T1>
                                       && (T1::Prec == Float32 || T2::Prec == Float64)
                                       && U1::isStrided()
                                       && U2::isStrided();
        };
    }
    /**
     * \class RValueMatrix: The base class of all matrixes
     */
    template<class Derived, Scalar ScalarT>
    class RValueMatrix : public CRTPBase<RValueMatrix<Derived, ScalarT>> {
        static_assert(!DeviceObj<Derived>, "[Error]: device_obj<> must be outside RValueMatrix<>");
        using This = RValueMatrix<Derived, ScalarT>;
        using Base = CRTPBase<This>;
    public:
        using ScalarType = ScalarT;
    protected:
        using T = ScalarType;
        using Tr = T::RealType;
        using Tv = T::ValueType;
        using Trv = Tr::ValueType;
        using Tc = T::ComplexType;
        using Tcv = Tc::ValueType;
    public:
        ~RValueMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        [[nodiscard]] auto operator*(this auto&&, Vector auto&& v) noexcept;
        [[nodiscard]] auto operator*(this auto&&, Matrix auto&& m) noexcept;
        [[nodiscard, gnu::always_inline]] auto operator-(this auto&&) noexcept;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Matrix auto&& target) const noexcept;
        template<ExecutePolicy P = Sequential>
        void assign_base(Matrix auto&& target) const noexcept;
        template<ExecutePolicy P = Sequential>
        void assign_add(Matrix auto&& target) const noexcept;
        void assert_assign(const Matrix auto& source) const noexcept;
        void assert_assign_lapack(const Matrix auto& source) const noexcept;

        [[nodiscard]] decltype(auto) calc(size_t row, size_t col) const;
        [[nodiscard]] decltype(auto) calc_value(size_t row, size_t col) const;
        [[nodiscard]] decltype(auto) calcFromMajorMinor(size_t major, size_t minor) const;
        void reverse(this const auto&, const Matrix auto& grad) noexcept;
        void reverse(this const auto&, const Matrix auto& y, const Matrix auto& grad) noexcept;

        void resize(const Matrix auto& m, auto&&... args);
        decltype(auto) resize(size_t r, size_t c, auto&&... args);

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
        [[nodiscard]] auto diag(this auto&&, ssize_t offset) noexcept;
        [[nodiscard]] auto triu(this auto&&) noexcept;
        [[nodiscard]] auto triu_unit(this auto&&) noexcept;
        [[nodiscard]] auto tril(this auto&&) noexcept;
        [[nodiscard]] auto tril_unit(this auto&&) noexcept;

        [[nodiscard]] Tr norm1() const;
        template<ExecutePolicy P = Sequential>
        [[nodiscard]] Tr norm1_power(unsigned int maxIteration) const;
        [[nodiscard]] CoDiff<Tr> normF() const;
        [[nodiscard]] Tr normInf() const;
        [[nodiscard]] T cond2() const;

        [[nodiscard]] Index2D argmax() const noexcept;
        [[nodiscard]] Index2D argmin() const noexcept;
        [[nodiscard]] T max() const;
        [[nodiscard]] T min() const;
        [[nodiscard]] CoDiff<T> sum() const;
        [[nodiscard]] auto sum_rows() const;
        [[nodiscard]] auto sum_cols() const;
        [[nodiscard]] CoDiff<T> mean() const;
        [[nodiscard]] T trace() const;
        [[nodiscard]] CoDiff<T> lnSumExp() const;
        [[nodiscard]] CoDiff<T> det() const;
        [[nodiscard]] Tr lnAbsDet() const;
        [[nodiscard]] auto sgndet() const;

        [[nodiscard]] auto format() const noexcept;
        [[nodiscard]] auto inv(this auto&&) noexcept;
        [[nodiscard]] auto pinv() const noexcept;
        [[nodiscard]] auto transpose(this auto&&) noexcept;
        [[nodiscard]] decltype(auto) conjugate(this auto&&) noexcept;
        [[nodiscard]] auto hermite(this auto&&) noexcept;
        [[nodiscard]] auto flatten(this auto&&) noexcept;

        [[nodiscard]] decltype(auto) reals(this auto&&) noexcept;
        [[nodiscard]] auto imags(this auto&&) noexcept;
        [[nodiscard]] auto squaredNorms(this auto&&) noexcept;
        [[nodiscard]] auto norms(this auto&&) noexcept;
        [[nodiscard]] decltype(auto) values(this auto&&) noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] decltype(auto) grads(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] auto getRow() const noexcept;
        [[nodiscard]] auto getCol() const noexcept;
        [[nodiscard]] size_t getSize() const noexcept { return getRow() * getCol(); }
        [[nodiscard]] auto getOrder() const noexcept;
        [[nodiscard]] size_t getMaxMajor(this const auto&) noexcept;
        [[nodiscard]] size_t getMaxMinor(this const auto&) noexcept;
        [[nodiscard]] bool empty() const noexcept { return Base::getDerived().getSize() == 0; }

        [[nodiscard]] constexpr bool isOverdetermined() const noexcept;
        [[nodiscard]] constexpr bool isUnderdetermined() const noexcept;
        [[nodiscard]] constexpr bool isSquare() const noexcept;
        [[nodiscard]] constexpr bool isSymm() const noexcept;
        [[nodiscard]] constexpr bool isHermite() const noexcept;
        [[nodiscard]] bool isFinite() const noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isComplex() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isDiffable() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isForwardDiff() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isReverseDiff() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isLValueMatrix() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isStrided() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isCompact() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isSparse() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isStaticSymm() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isStaticHermite() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isStaticSquare() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isColMatrix() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isRowMatrix() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isBothMajor() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isFastAssign() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getRowAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getColAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static auto getMajor() noexcept;
        [[nodiscard]] static size_t rowFromMajorMinor(size_t major, size_t minor) noexcept { return MatrixMajor::rowFromMajorMinor<Derived>(major, minor); }
        [[nodiscard]] static size_t colFromMajorMinor(size_t major, size_t minor) noexcept { return MatrixMajor::colFromMajorMinor<Derived>(major, minor); }
        __host__ __device__ consteval static void static_assert_assign(const Scalar auto& source) noexcept;
        __host__ __device__ consteval static void static_assert_assign(const Matrix auto& source) noexcept;
    protected:
        RValueMatrix() = default;
        RValueMatrix(const This&) = default;
        RValueMatrix(This&&) noexcept = default;
        /* Static members */
        [[nodiscard]] consteval static int calcBlockingSize(int CacheSize) noexcept;
        __host__ __device__ static void checkBlock(const Matrix auto& m, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept;
    };

    template<Matrix M1, Matrix M2>
    bool matrixNear(const M1& m1, const M2& m2, double precision) noexcept {
        using T = Internal::BinaryScalarOpRtnTy<typename M1::ScalarType, typename M2::ScalarType>::Type;
        assert(m1.getRow() == m2.getRow());
        assert(m1.getCol() == m2.getCol());
        for (size_t i = 0; i < m1.getCol(); ++i)
            for (size_t j = 0; j < m1.getRow(); ++j)
                if (!scalarNear(T(m1.calc(j, i)), T(m2.calc(j, i)), precision))
                    return false;
        return true;
    }

    template<Matrix M1, Matrix M2>
    bool matrixNear(const M1& m1, const M2& m2, uint64_t ulp) noexcept {
        using T = Internal::BinaryScalarOpRtnTy<typename M1::ScalarType, typename M2::ScalarType>::Type;
        assert(m1.getRow() == m2.getRow());
        assert(m1.getCol() == m2.getCol());
        for (size_t i = 0; i < m1.getCol(); ++i)
            for (size_t j = 0; j < m1.getRow(); ++j)
                if (!scalarNear(T(m1.calc(j, i)), T(m2.calc(j, i)), ulp))
                    return false;
        return true;
    }

    bool operator==(const Matrix auto& m1, const Matrix auto& m2) noexcept {
        if (m1.getRow() != m2.getRow())
            return false;
        if (m1.getCol() != m2.getCol())
            return false;
        for (size_t major = 0; major < m1.getMaxMajor(); ++major) {
            for (size_t minor = 0; minor < m1.getMaxMinor(); ++minor) {
                size_t r = m1.rowFromMajorMinor(major, minor);
                size_t c = m1.colFromMajorMinor(major, minor);
                if (m1.calc(r, c) != m2.calc(r, c))
                    return false;
            }
        }
        return true;
    }

    bool operator!=(const Matrix auto& m1, const Matrix auto& m2) noexcept { return !(m1 == m2); }

    std::ostream& operator<<(std::ostream& os, const Matrix auto& m) noexcept {
        return os << std::format("{}", m.format());
    }
}

namespace Physica {
    template<class T, Scalar S>
    class Traits<RValueMatrix<T, S>> {
    public:
        using Derived = T;
    };
}

#include "RValueMatrixImpl/RValueMatrixImpl.h"
#include "RValueMatrixImpl/MainDiag.h"
#include "RValueMatrixImpl/OffsetDiag.h"
#include "MatrixProduct/GEMM.h"
#include "MatrixProduct/GEMV.h"
#include "MatrixProduct/GEVM.h"
#include "MatrixProduct/Kronecker.h"
#include "RValueMatrixImpl/MatrixNorm.h"
#include "MatrixExpr.h"
