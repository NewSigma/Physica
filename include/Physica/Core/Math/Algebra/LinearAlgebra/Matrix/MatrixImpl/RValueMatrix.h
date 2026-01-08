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
    template<class Derived> class ContinuousMatrix;
    template<class, bool ReduceCol> class MatrixSum;
    template<Matrix> class DiagVectorR;
    template<class> class Inverse;
    template<Matrix> class PseudoInverse;
    template<class> class Transpose;
    template<class> class Conjugate;
    template<class> class Hermite;
    template<class> class FlattenR;
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
        class EnableMKL<M1, M2> {
            using U1 = std::remove_cvref<M1>::type;
            using U2 = std::remove_cvref<M2>::type;
            using ScalarType1 = U1::ScalarType;
            using ScalarType2 = U2::ScalarType;
        public:
            constexpr static bool value = HasMKL()
                                       && std::same_as<ScalarType1, ScalarType2>
                                       && (ScalarType1::Prec == Float32 || ScalarType1::Prec == Float64)
                                       && is_continuous<U1>::value
                                       && is_continuous<U2>::value
                                       && !Diffable<U1>;
        };
    }
    /**
     * \class RValueMatrix: The base class of all matrixes
     */
    template<class Derived>
    class RValueMatrix : public CRTPBase<RValueMatrix<Derived>> {
        static_assert(!CUDA<Derived>, "[Error]: device_obj<> must be outside RValueMatrix<>");
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

        using Tm = std::conditional<isComplex, typename Tc::MKL_Complex, typename T::MachineType>::type;
    public:
        ~RValueMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        template<Vector V>
        [[nodiscard]] auto operator*(V&& v) const& noexcept requires(RowAtCompile != 1 && !CUDA<V>);
        template<Vector V>
        [[nodiscard]] auto operator*(V&& v) && noexcept requires(RowAtCompile != 1 && !CUDA<V>);
        template<Vector V>
        [[nodiscard]] auto operator*(const V& v) const noexcept requires(RowAtCompile == 1 && !CUDA<V>);
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Matrix auto&& target) const;
        template<ExecutePolicy P = Sequential>
        void assign_base(Matrix auto&& target) const;
        void assign_add(Matrix auto&& target) const;
        void assert_assign(const Matrix auto& source) const noexcept;
        void assert_assign_mkl(const Matrix auto& source) const noexcept;

        [[nodiscard]] decltype(auto) calc(size_t row, size_t col) const { return Base::getDerived().calc(row, col); }
        [[nodiscard]] decltype(auto) calc_value(size_t row, size_t col) const { return Base::getDerived().calc_value(row, col); }
        [[nodiscard]] decltype(auto) calcFromMajorMinor(size_t major, size_t minor) const;
        void reverse(const Matrix auto& y, const Matrix auto& grad) const noexcept;

        void resize(const Matrix auto& m, auto&&... args);
        decltype(auto) resize(size_t r, size_t c, auto&&... args);

        [[nodiscard]] auto row(size_t r) noexcept;
        [[nodiscard]] const auto row(size_t r) const noexcept;
        [[nodiscard]] auto col(size_t c) noexcept;
        [[nodiscard]] const auto col(size_t c) const noexcept;
        [[nodiscard]] auto rows(size_t fromRow, size_t rowCount) noexcept;
        [[nodiscard]] const auto rows(size_t fromRow, size_t rowCount) const noexcept;
        [[nodiscard]] auto topRows(size_t to) noexcept;
        [[nodiscard]] const auto topRows(size_t to) const noexcept;
        [[nodiscard]] auto bottomRows(size_t from) noexcept;
        [[nodiscard]] const auto bottomRows(size_t from) const noexcept;
        [[nodiscard]] auto cols(size_t fromCol, size_t colCount) noexcept;
        [[nodiscard]] const auto cols(size_t fromCol, size_t colCount) const noexcept;
        [[nodiscard]] auto leftCols(size_t to) noexcept;
        [[nodiscard]] const auto leftCols(size_t to) const noexcept;
        [[nodiscard]] auto rightCols(size_t from) noexcept;
        [[nodiscard]] const auto rightCols(size_t from) const noexcept;
        [[nodiscard]] auto topLeftCorner(size_t toRow, size_t toCol) noexcept;
        [[nodiscard]] const auto topLeftCorner(size_t toRow, size_t toCol) const noexcept;
        [[nodiscard]] auto topLeftCorner(size_t to) noexcept;
        [[nodiscard]] const auto topLeftCorner(size_t to) const noexcept;
        [[nodiscard]] auto topRightCorner(size_t toRow, size_t fromCol) noexcept;
        [[nodiscard]] const auto topRightCorner(size_t toRow, size_t fromCol) const noexcept;
        [[nodiscard]] auto bottomLeftCorner(size_t fromRow, size_t toCol) noexcept;
        [[nodiscard]] const auto bottomLeftCorner(size_t fromRow, size_t toCol) const noexcept;
        [[nodiscard]] auto bottomRightCorner(size_t fromRow, size_t fromCol) noexcept;
        [[nodiscard]] const auto bottomRightCorner(size_t fromRow, size_t fromCol) const noexcept;
        [[nodiscard]] auto bottomRightCorner(size_t from) noexcept;
        [[nodiscard]] const auto bottomRightCorner(size_t from) const noexcept;
        [[nodiscard]] auto block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept;
        [[nodiscard]] const auto block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) const noexcept;
        [[nodiscard]] auto diag(this auto&&) noexcept;
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
        [[nodiscard]] auto inv() const noexcept;
        [[nodiscard]] auto pinv() const noexcept;
        [[nodiscard]] auto transpose() const noexcept;
        [[nodiscard]] decltype(auto) conjugate() const noexcept;
        [[nodiscard]] auto hermite() const noexcept;
        [[nodiscard]] auto flatten() const noexcept;

        [[nodiscard]] decltype(auto) reals() const noexcept;
        [[nodiscard]] auto imags() const noexcept;
        [[nodiscard]] auto squaredNorms() const noexcept;
        [[nodiscard]] auto norms() const noexcept;
        [[nodiscard]] decltype(auto) values(this auto&&) noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] decltype(auto) grads() const noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return Base::getDerived().getRow(); }
        [[nodiscard]] size_t getCol() const noexcept { return Base::getDerived().getCol(); }
        [[nodiscard]] size_t getSize() const noexcept { return getRow() * getCol(); }
        [[nodiscard]] size_t getMaxMajor() const noexcept { return MatrixOption::getMaxMajor<Derived>(Base::getDerived()); }
        [[nodiscard]] size_t getMaxMinor() const noexcept { return MatrixOption::getMaxMinor<Derived>(Base::getDerived()); }

        [[nodiscard]] bool isOverdetermined() const noexcept;
        [[nodiscard]] bool isUnderdetermined() const noexcept;
        [[nodiscard]] bool isSquare() const noexcept;
        [[nodiscard]] bool isSymm() const noexcept;
        [[nodiscard]] bool isHermite() const noexcept;
        [[nodiscard]] bool isFinite() const noexcept;
        /* Static members */
        [[nodiscard]] static size_t rowFromMajorMinor(size_t major, size_t minor) noexcept { return MatrixOption::rowFromMajorMinor<Derived>(major, minor); }
        [[nodiscard]] static size_t colFromMajorMinor(size_t major, size_t minor) noexcept { return MatrixOption::colFromMajorMinor<Derived>(major, minor); }
        __host__ __device__ static void static_assert_assign(const Matrix auto& source) noexcept;
    protected:
        RValueMatrix() = default;
        RValueMatrix(const This&) = default;
        RValueMatrix(This&&) noexcept = default;
        /* Operations */
        template<int GradOrder>
        auto grads_impl() const noexcept;
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

    template<Matrix M>
    bool operator==(const M& m1, const M& m2) {
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

    template<Matrix M>
    bool operator!=(const M& m1, const M& m2) { return !(m1 == m2); }

    std::ostream& operator<<(std::ostream& os, const Matrix auto& m) noexcept {
        return os << std::format("{}", m.format());
    }
}

namespace Physica {
    template<class T>
    class Traits<RValueMatrix<T>> : public Traits<T> {
    public:
        using Derived = T;
    };
}

#include "RValueMatrixImpl/RValueMatrixImpl.h"
#include "RValueMatrixImpl/DiagVector.h"
#include "MatrixProduct/GEMM.h"
#include "MatrixProduct/GEMV.h"
#include "MatrixProduct/GEVM.h"
#include "MatrixProduct/Kronecker.h"
#include "RValueMatrixImpl/MatrixNorm.h"
#include "MatrixExpr.h"
