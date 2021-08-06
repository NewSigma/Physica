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

    template<class MatrixType1, class MatrixType2>
    template<class OtherDerived>
    void MatrixProduct<MatrixType1, MatrixType2>::assignTo(LValueMatrix<OtherDerived>& target) const {
        constexpr static int defaultMajor = Internal::ProductOption<MatrixType1, MatrixType2>::Major;
        constexpr static bool isDynamic = defaultMajor == Dynamic;
        using TargetType = LValueMatrix<OtherDerived>;

        if constexpr (isDynamic) {
            for (size_t i = 0; i < target.getMaxMajor(); ++i)
                for (size_t j = 0; j < target.getMaxMinor(); ++j)
                    target.getElementFromMajorMinor(i, j) = calc(TargetType::rowFromMajorMinor(i, j),
                                                                 TargetType::columnFromMajorMinor(i, j));
        }
        else {
            for (size_t i = 0; i < (defaultMajor == DenseMatrixOption::Column ? getColumn() : getRow()); ++i) {
                for (size_t j = 0; j < (defaultMajor == DenseMatrixOption::Column ?  getRow() : getColumn()); ++j) {
                    const size_t r = DefaultType::rowFromMajorMinor(i, j);
                    const size_t c = DefaultType::columnFromMajorMinor(i, j);
                    target(r, c) = calc(r, c);
                }
            }
        }
    }

    template<class T1, class T2>
    typename MatrixProduct<T1, T2>::ScalarType MatrixProduct<T1, T2>::calc(size_t row, size_t column) const {
        ScalarType result = 0;
        for (size_t i = 0; i < mat1.getColumn(); ++i)
            result += ScalarType(mat1(row, i) * mat2(i, column));
        return result;
    }

    template<class Derived, class OtherDerived>
    inline MatrixProduct<Derived, OtherDerived>
    operator*(const RValueMatrix<Derived>& mat1, const RValueMatrix<OtherDerived>& mat2) {
        assert(mat1.getColumn() == mat2.getRow());
        return MatrixProduct(mat1, mat2);
    }
}
