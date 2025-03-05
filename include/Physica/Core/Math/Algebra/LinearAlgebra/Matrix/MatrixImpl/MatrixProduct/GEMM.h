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
    template<Scalar T,
             int Option,
             size_t Row,
             size_t Col,
             class Allocator>
    class DenseMatrix;

    namespace Internal {
        template<Matrix T1, Matrix T2>
        struct ProductOption {
            constexpr static bool SameMajor = MatrixOption::isSameMajor<T1, T2>();
            constexpr static bool RowMajor = MatrixOption::isRowMatrix<T1>();
            constexpr static int Major = SameMajor ? (RowMajor ? int(MatrixOption::Col)
                                                               : int(MatrixOption::Row))
                                                   : int(MatrixOption::AnyMajor);
            constexpr static int Storage = (MatrixOption::isElementMatrix<T1>() && MatrixOption::isElementMatrix<T2>())
                                         ? MatrixOption::Element
                                         : MatrixOption::Vector;
            constexpr static int Option = (Major == MatrixOption::AnyMajor ? MatrixOption::Col : Major) | Storage;
        };
    }

    template<Matrix T1, Matrix T2>
    class MatrixProduct : public RValueMatrix<MatrixProduct<T1, T2>> {
        using This = MatrixProduct<T1, T2>;
        using Base = RValueMatrix<This>;
    public:
        using typename Base::ScalarType;
        using Base::isReverseDiff;
        using DefaultType = DenseMatrix<ScalarType,
                                        Internal::ProductOption<T1, T2>::Option,
                                        Base::RowAtCompile,
                                        Base::ColAtCompile,
                                        HostAllocator<ScalarType>>;
    private:
        const T1& mat1;
        const T2& mat2;
    public:
        MatrixProduct(const T1& mat1_, const T2& mat2_) : mat1(mat1_), mat2(mat2_) {
            assert(mat1.getCol() == mat2.getRow());
        }
        MatrixProduct(const This&) = default;
        MatrixProduct(This&&) noexcept = default;
        ~MatrixProduct() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<Matrix M>
        void assign(LValueMatrix<M>& target) const;
        [[nodiscard]] DefaultType compute() const { return DefaultType(*this); }

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const;

        template<Matrix M>
        void reverse(const M& grad) const noexcept requires(isReverseDiff);
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat1.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const { return mat2.getCol(); }
        [[nodiscard]] const T1& getLHS() const noexcept { return mat1; }
        [[nodiscard]] const T2& getRHS() const noexcept { return mat2; }
    };

    template<Matrix T1, Matrix T2>
    template<Matrix M>
    void MatrixProduct<T1, T2>::assign(LValueMatrix<M>& target) const {
        constexpr static int defaultMajor = Internal::ProductOption<T1, T2>::Major;
        constexpr static bool isAnyMajor = defaultMajor == MatrixOption::AnyMajor;
        using TargetType = LValueMatrix<M>;
        using T = TargetType::ScalarType;
        if constexpr (isAnyMajor) {
            for (size_t i = 0; i < target.getMaxMajor(); ++i) {
                for (size_t j = 0; j < target.getMaxMinor(); ++j) {
                    const size_t r = MatrixOption::rowFromMajorMinor<TargetType>(i, j);
                    const size_t c = MatrixOption::colFromMajorMinor<TargetType>(i, j);
                    target.refFromMajorMinor(i, j) = T(calc(r, c));
                }
            }
        }
        else {
            for (size_t i = 0; i < (defaultMajor == MatrixOption::Col ? getCol() : getRow()); ++i) {
                for (size_t j = 0; j < (defaultMajor == MatrixOption::Col ?  getRow() : getCol()); ++j) {
                    const size_t r = MatrixOption::rowFromMajorMinor<DefaultType>(i, j);
                    const size_t c = MatrixOption::colFromMajorMinor<DefaultType>(i, j);
                    target(r, c) = T(calc(r, c));
                }
            }
        }
    }

    template<Matrix T1, Matrix T2>
    auto MatrixProduct<T1, T2>::calc(size_t row, size_t col) const -> ScalarType {
        ScalarType result(0);
        for (size_t i = 0; i < mat1.getCol(); ++i)
            result += ScalarType(mat1.calc(row, i)) * ScalarType(mat2.calc(i, col));
        return result;
    }

    template<Matrix T1, Matrix T2>
    template<Matrix M>
    void MatrixProduct<T1, T2>::reverse(const M& grad) const noexcept requires(isReverseDiff) {
        if constexpr (ReverseDiff<T1>)
            mat1.reverse(grad * mat2.transpose());
        if constexpr (ReverseDiff<T2>)
            mat2.reverse(mat1.transpose() * grad);
    }

    template<Matrix T1, Matrix T2>
    [[nodiscard]] inline auto operator*(const T1& mat1, const T2& mat2) noexcept
            requires(((T1::ColAtCompile != 1 && T2::ColAtCompile != 1) || (T1::ColAtCompile == 1 && T2::ColAtCompile == 1)) && !CUDA<T1> && !CUDA<T2>) {
        assert(mat1.getCol() == mat2.getRow());
        return MatrixProduct<T1, T2>(mat1, mat2);
    }
}

namespace Physica {
    template<Matrix T1, Matrix T2>
    class Traits<MatrixProduct<T1, T2>> {
        static_assert(T1::ColAtCompile == T2::RowAtCompile ||
                      T1::ColAtCompile == Dynamic ||
                      T2::RowAtCompile == Dynamic,
                      "[Error]: Row and column do not match in matrix-vector product");
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename T1::ScalarType, typename T2::ScalarType>::Type;
        constexpr static int Option = MatrixOption::AnyMajor | MatrixOption::AnyStorage;
        constexpr static size_t RowAtCompile = T1::RowAtCompile;
        constexpr static size_t ColAtCompile = T2::ColAtCompile;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;
    };
}
