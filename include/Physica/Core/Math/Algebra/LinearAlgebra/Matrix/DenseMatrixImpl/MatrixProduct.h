/*
 * Copyright 2021 WeiBo He.
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

namespace Physica::Core {
    template<class MatrixType1, class MatrixType2> class MatrixProduct;
    template<class VectorType, class MatrixType> class VectorMatrixProduct;
    template<class MatrixType, class VectorType> class MatrixVectorProduct;

    namespace Internal {
        template<class MatrixType1, class MatrixType2>
        struct ProductOption {
            constexpr static bool SameMajor = DenseMatrixOption::isSameMajor<MatrixType1, MatrixType2>();
            constexpr static bool RowMajor = DenseMatrixOption::isRowMatrix<MatrixType1>();
            constexpr static int Major = SameMajor ? (RowMajor ? int(DenseMatrixOption::Column)
                                                               : int(DenseMatrixOption::Row))
                                                   : Dynamic;
            constexpr static int Storage = (DenseMatrixOption::isElementMatrix<MatrixType1>() && DenseMatrixOption::isElementMatrix<MatrixType2>())
                                         ? DenseMatrixOption::Element
                                         : DenseMatrixOption::Vector;
            constexpr static int MatrixOption = (Major == Dynamic ? DenseMatrixOption::Column : Major) | Storage;
        };

        template<class MatrixType1, class MatrixType2>
        class Traits<MatrixProduct<MatrixType1, MatrixType2>> {
        public:
            using ScalarType = typename Internal::BinaryScalarOpReturnType<typename MatrixType1::ScalarType,
                                                                           typename MatrixType2::ScalarType>::Type;
            constexpr static size_t RowAtCompile = MatrixType1::RowAtCompile;
            constexpr static size_t ColumnAtCompile = MatrixType2::ColumnAtCompile;
            constexpr static size_t MaxRowAtCompile = MatrixType1::MaxRowAtCompile;
            constexpr static size_t MaxColumnAtCompile = MatrixType2::MaxColumnAtCompile;
            constexpr static size_t SizeAtCompile = RowAtCompile * ColumnAtCompile;
            constexpr static size_t MaxSizeAtCompile = SizeAtCompile;
        };

        template<class VectorType, class MatrixType>
        class Traits<VectorMatrixProduct<VectorType, MatrixType>> {
        public:
            using ScalarType = typename Internal::BinaryScalarOpReturnType<typename VectorType::ScalarType,
                                                                           typename MatrixType::ScalarType>::Type;
            constexpr static size_t RowAtCompile = VectorType::SizeAtCompile;
            constexpr static size_t ColumnAtCompile = MatrixType::ColumnAtCompile;
            constexpr static size_t MaxRowAtCompile = VectorType::MaxSizeAtCompile;
            constexpr static size_t MaxColumnAtCompile = MatrixType::MaxColumnAtCompile;
            constexpr static size_t SizeAtCompile = RowAtCompile * ColumnAtCompile;
            constexpr static size_t MaxSizeAtCompile = SizeAtCompile;
        };

        template<class MatrixType, class VectorType>
        class Traits<MatrixVectorProduct<MatrixType, VectorType>> {
        public:
            using ScalarType = typename Internal::BinaryScalarOpReturnType<typename MatrixType::ScalarType,
                                                                           typename VectorType::ScalarType>::Type;
            constexpr static size_t SizeAtCompile = MatrixType::RowAtCompile;
            constexpr static size_t MaxSizeAtCompile = MatrixType::MaxRowAtCompile;
            using PacketType = typename Internal::BestPacket<ScalarType, SizeAtCompile>::Type;
        };
    }

    template<class MatrixType1, class MatrixType2>
    class MatrixProduct : public RValueMatrix<MatrixProduct<MatrixType1, MatrixType2>> {
        static_assert(MatrixType1::ColumnAtCompile == MatrixType2::RowAtCompile ||
                      MatrixType1::ColumnAtCompile == Dynamic ||
                      MatrixType2::RowAtCompile == Dynamic,
                      "Row and column do not match in matrix product");
    public:
        using Base = RValueMatrix<MatrixProduct<MatrixType1, MatrixType2>>;
        using typename Base::ScalarType;
        using DefaultType = DenseMatrix<ScalarType,
                                        Internal::ProductOption<MatrixType1, MatrixType2>::MatrixOption,
                                        Base::RowAtCompile,
                                        Base::ColumnAtCompile,
                                        Base::MaxRowAtCompile,
                                        Base::MaxColumnAtCompile>;
    private:
        const MatrixType1& mat1;
        const MatrixType2& mat2;
    public:
        MatrixProduct(const RValueMatrix<MatrixType1>& mat1_, const RValueMatrix<MatrixType2>& mat2_)
                : mat1(mat1_.getDerived()), mat2(mat2_.getDerived()) {
            assert(mat1.getColumn() == mat2.getRow());
        }
        /* Operations */
        template<class OtherDerived>
        void assignTo(LValueMatrix<OtherDerived>& target) const;
        [[nodiscard]] DefaultType compute() const { return DefaultType(*this); }
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t column) const;
        [[nodiscard]] size_t getRow() const { return mat1.getRow(); }
        [[nodiscard]] size_t getColumn() const { return mat2.getColumn(); }
        [[nodiscard]] const MatrixType1& getLHS() const noexcept { return mat1; }
        [[nodiscard]] const MatrixType2& getRHS() const noexcept { return mat2; }
    };

    template<class VectorType, class MatrixType>
    class VectorMatrixProduct : public RValueMatrix<VectorMatrixProduct<VectorType, MatrixType>> {
        static_assert(MatrixType::RowAtCompile == 1 || MatrixType::RowAtCompile == Dynamic,
                      "Row and column do not match in matrix product");
    public:
        using Base = RValueMatrix<VectorMatrixProduct<VectorType, MatrixType>>;
        using typename Base::ScalarType;
    private:
        const VectorType& vec;
        const MatrixType& mat;
    public:
        VectorMatrixProduct(const RValueVector<VectorType>& vec_, const RValueMatrix<MatrixType>& mat_)
                : vec(vec_.getDerived()), mat(mat_.getDerived()) {
            assert(mat.getRow() == 1);
        }
        /* Operations */
        template<class OtherDerived>
        void assignTo(LValueMatrix<OtherDerived>& target) const;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t column) const;
        [[nodiscard]] size_t getRow() const { return vec.getLength(); }
        [[nodiscard]] size_t getColumn() const { return mat.getColumn(); }
        [[nodiscard]] const VectorType& getLHS() const noexcept { return vec; }
        [[nodiscard]] const MatrixType& getRHS() const noexcept { return mat; }
    };

    template<class MatrixType, class VectorType>
    class MatrixVectorProduct : public RValueVector<MatrixVectorProduct<MatrixType, VectorType>> {
        static_assert(MatrixType::ColumnAtCompile == VectorType::SizeAtCompile ||
                      MatrixType::ColumnAtCompile == Dynamic ||
                      VectorType::SizeAtCompile == Dynamic,
                      "Row and column do not match in matrix product");
    public:
        using Base = RValueVector<MatrixVectorProduct<MatrixType, VectorType>>;
        using typename Base::ScalarType;
    private:
        const MatrixType& mat;
        const VectorType& vec;
    public:
        MatrixVectorProduct(const RValueMatrix<MatrixType>& mat_, const RValueVector<VectorType>& vec_)
                : mat(mat_.getDerived()), vec(vec_.getDerived()) {
            assert(mat.getColumn() == vec.getLength());
        }
        /* Operations */
        template<class OtherDerived>
        void assignTo(LValueVector<OtherDerived>& target) const;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t index) const;
        [[nodiscard]] size_t getLength() const { return mat.getRow(); }
        [[nodiscard]] const MatrixType& getLHS() const noexcept { return mat; }
        [[nodiscard]] const VectorType& getRHS() const noexcept { return vec; }
    };

    template<class MatrixType1, class MatrixType2>
    inline typename std::enable_if<MatrixType2::ColumnAtCompile != 1, MatrixProduct<MatrixType1, MatrixType2>>::type
    operator*(const RValueMatrix<MatrixType1>& mat1, const RValueMatrix<MatrixType2>& mat2) {
        assert(mat1.getColumn() == mat2.getRow());
        return MatrixProduct(mat1, mat2);
    }

    template<class VectorType, class MatrixType>
    inline VectorMatrixProduct<VectorType, MatrixType>
    operator*(const RValueVector<VectorType>& vec, const RValueMatrix<MatrixType>& mat) {
        assert(mat.getRow() == 1);
        return VectorMatrixProduct(vec, mat);
    }

    template<class MatrixType, class VectorType>
    inline MatrixVectorProduct<MatrixType, VectorType>
    operator*(const RValueMatrix<MatrixType>& mat, const RValueVector<VectorType>& vec) {
        assert(mat.getColumn() == vec.getLength());
        return MatrixVectorProduct(mat, vec);
    }
}

#include "MatrixProductImpl.h"
