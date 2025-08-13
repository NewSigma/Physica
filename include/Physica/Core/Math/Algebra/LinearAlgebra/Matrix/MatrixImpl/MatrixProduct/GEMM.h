/*
 * Copyright 2024-2025 Weibo He.
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

#include "../RValueMatrix.h"

namespace Physica {
    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    class DenseMatrix;

    namespace Internal {
        template<Matrix M1, Matrix M2>
        struct ProductOption {
            constexpr static bool SameMajor = MatrixOption::isSameMajor<M1, M2>();
            constexpr static bool RowMajor = MatrixOption::isRowMatrix<M1>();
            constexpr static int Major = SameMajor ? (RowMajor ? int(MatrixOption::Col)
                                                               : int(MatrixOption::Row))
                                                   : int(MatrixOption::AnyMajor);
            constexpr static int Storage = (MatrixOption::isElementMatrix<M1>() && MatrixOption::isElementMatrix<M2>())
                                         ? MatrixOption::Element
                                         : MatrixOption::Vector;
            constexpr static int Option = (Major == MatrixOption::AnyMajor ? MatrixOption::Col : Major) | Storage;
        };
    }

    template<Matrix M1, Matrix M2>
    class GEMM : public RValueMatrix<GEMM<M1, M2>> {
        using This = GEMM<M1, M2>;
        using Base = RValueMatrix<This>;
    public:
        using Base::isReverseDiff;
        using Base::isComplex;
    protected:
        using typename Base::T;
    private:
        using DefaultType = DenseMatrix<T,
                                        Internal::ProductOption<M1, M2>::Option,
                                        Base::RowAtCompile,
                                        Base::ColAtCompile,
                                        HostAllocator<T>>;

        const M1& mat1;
        const M2& mat2;
    public:
        GEMM(const M1& mat1_, const M2& mat2_);
        GEMM(const This&) = default;
        GEMM(This&&) noexcept = default;
        ~GEMM() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(Matrix auto& target) const;
        void assign_base(Matrix auto& target) const;
        void assign_mkl(Matrix auto& target) const noexcept;
        [[nodiscard]] DefaultType compute() const { return DefaultType(*this); }

        [[nodiscard]] T calc(size_t row, size_t col) const;

        void reverse(const Matrix auto& grad) const noexcept requires(isReverseDiff);
        /* Getters */
        [[nodiscard]] size_t getRow() const { return mat1.getRow(); }
        [[nodiscard]] size_t getCol() const { return mat2.getCol(); }
        [[nodiscard]] const M1& getLHS() const noexcept { return mat1; }
        [[nodiscard]] const M2& getRHS() const noexcept { return mat2; }
        /* Friends */
        friend class device_obj<This>;
    };

    template<Matrix M1, Matrix M2>
    GEMM<M1, M2>::GEMM(const M1& mat1_, const M2& mat2_) : mat1(mat1_), mat2(mat2_) {
        assert(mat1.getRow() > 0);
        assert(mat1.getCol() > 0);
        assert(mat2.getCol() > 0);
        assert(mat1.getCol() == mat2.getRow());
    }

    template<Matrix M1, Matrix M2>
    void GEMM<M1, M2>::assign(Matrix auto& target) const {
        using M = std::remove_cvref_t<decltype(target)>;
        constexpr bool GoodScalar = T::Prec == Float32 || T::Prec == Float64;
        constexpr bool SameScalar = std::same_as<typename M1::ScalarType, typename M2::ScalarType> && std::same_as<T, typename M::ScalarType>;
        constexpr bool SameMajor = MatrixOption::getMajor<M1>() == MatrixOption::getMajor<M2>();
        constexpr bool isContinuous = is_continuous<M1>::value && is_continuous<M2>::value && is_continuous<M>::value;
        constexpr bool isElement = MatrixOption::isElementMatrix<M1>() && MatrixOption::isElementMatrix<M2>() && MatrixOption::isElementMatrix<M>();
        constexpr bool isDiffable = Diffable<T>;
        constexpr bool UseMKL = HasMKL() && GoodScalar && SameScalar && SameMajor && isContinuous && isElement && !isDiffable;
        constexpr bool SmallMatrix1 = 0 < M1::SizeAtCompile && M1::SizeAtCompile <= 64;
        constexpr bool SmallMatrix2 = 0 < M2::SizeAtCompile && M2::SizeAtCompile <= 64;
        constexpr bool SmallMatrix = SmallMatrix1 && SmallMatrix2;
        assert(target.getRow() == getRow() && target.getCol() == getCol());
        if constexpr (UseMKL && !SmallMatrix)
            assign_mkl(target);
        else
            assign_base(target);
    }

    template<Matrix M1, Matrix M2>
    void GEMM<M1, M2>::assign_base(Matrix auto& target) const {
        constexpr static int defaultMajor = Internal::ProductOption<M1, M2>::Major;
        constexpr static bool isAnyMajor = defaultMajor == MatrixOption::AnyMajor;
        if constexpr (isAnyMajor) {
            for (size_t i = 0; i < target.getMaxMajor(); ++i) {
                for (size_t j = 0; j < target.getMaxMinor(); ++j) {
                    const size_t r = target.rowFromMajorMinor(i, j);
                    const size_t c = target.colFromMajorMinor(i, j);
                    target.refFromMajorMinor(i, j) = calc(r, c);
                }
            }
        }
        else {
            for (size_t i = 0; i < (defaultMajor == MatrixOption::Col ? getCol() : getRow()); ++i) {
                for (size_t j = 0; j < (defaultMajor == MatrixOption::Col ?  getRow() : getCol()); ++j) {
                    const size_t r = MatrixOption::rowFromMajorMinor<DefaultType>(i, j);
                    const size_t c = MatrixOption::colFromMajorMinor<DefaultType>(i, j);
                    target(r, c) = calc(r, c);
                }
            }
        }
    }

    template<Matrix M1, Matrix M2>
    auto GEMM<M1, M2>::calc(size_t row, size_t col) const -> T {
        T result(0);
        for (size_t i = 0; i < mat1.getCol(); ++i)
            result += mat1.calc(row, i) * mat2.calc(i, col);
        return result;
    }

    template<Matrix M1, Matrix M2>
    void GEMM<M1, M2>::reverse(const Matrix auto& grad) const noexcept requires(isReverseDiff) {
        if constexpr (ReverseDiff<M1>)
            mat1.reverse(grad * mat2.transpose());
        if constexpr (ReverseDiff<M2>)
            mat2.reverse(mat1.transpose() * grad);
    }
}

namespace Physica {
    template<Matrix M1, Matrix M2>
    class Traits<GEMM<M1, M2>> {
        static_assert(M1::ColAtCompile == M2::RowAtCompile ||
                      M1::ColAtCompile == Dynamic ||
                      M2::RowAtCompile == Dynamic,
                      "[Error]: Row and column do not match in matrix-vector product");
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename M1::ScalarType, typename M2::ScalarType>::Type;
        constexpr static int Option = MatrixOption::AnyMajor | MatrixOption::AnyStorage;
        constexpr static size_t RowAtCompile = M1::RowAtCompile;
        constexpr static size_t ColAtCompile = M2::ColAtCompile;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;
    };
}

#ifdef PHYSICA_MKL
    #include "GEMM_MKL.h"
#endif
