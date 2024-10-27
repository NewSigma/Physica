/*
 * Copyright 2024 Weibo He.
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
    namespace Internal {
        template<class MatrixType1, class MatrixType2>
        struct ProductOption {
            constexpr static bool SameMajor = MatrixOption::isSameMajor<MatrixType1, MatrixType2>();
            constexpr static bool RowMajor = MatrixOption::isRowMatrix<MatrixType1>();
            constexpr static int Major = SameMajor ? (RowMajor ? int(MatrixOption::Column)
                                                               : int(MatrixOption::Row))
                                                   : int(MatrixOption::AnyMajor);
            constexpr static int Storage = (MatrixOption::isElementMatrix<MatrixType1>() && MatrixOption::isElementMatrix<MatrixType2>())
                                         ? MatrixOption::Element
                                         : MatrixOption::Vector;
            constexpr static int Option = (Major == MatrixOption::AnyMajor ? MatrixOption::Column : Major) | Storage;
        };
    }

    template<class MatrixType1, class MatrixType2>
    class MatrixProduct : public RValueMatrix<MatrixProduct<MatrixType1, MatrixType2>> {
        using This = MatrixProduct<MatrixType1, MatrixType2>;
        using Base = RValueMatrix<This>;
    public:
        using Base::isReverseDiff;
        using typename Base::ScalarType;
        using DefaultType = DenseMatrix<ScalarType,
                                        Internal::ProductOption<MatrixType1, MatrixType2>::Option,
                                        Base::RowAtCompile,
                                        Base::ColumnAtCompile,
                                        HostAllocator<ScalarType>>;
    private:
        const MatrixType1& mat1;
        const MatrixType2& mat2;
    public:
        MatrixProduct(const RValueMatrix<MatrixType1>& mat1_, const RValueMatrix<MatrixType2>& mat2_)
                : mat1(mat1_.getDerived()), mat2(mat2_.getDerived()) {
            assert(mat1.getColumn() == mat2.getRow());
        }
        MatrixProduct(const This&) = delete;
        MatrixProduct(This&&) noexcept = delete;
        ~MatrixProduct() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<class OtherDerived>
        void assignTo(LValueMatrix<OtherDerived>& target) const;
        [[nodiscard]] DefaultType compute() const { return DefaultType(*this); }
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t column) const;
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat1.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return mat2.getColumn(); }
        [[nodiscard]] const MatrixType1& getLHS() const noexcept { return mat1; }
        [[nodiscard]] const MatrixType2& getRHS() const noexcept { return mat2; }
    };

    template<class MatrixType1, class MatrixType2>
    template<class OtherDerived>
    void MatrixProduct<MatrixType1, MatrixType2>::assignTo(LValueMatrix<OtherDerived>& target) const {
        constexpr static int defaultMajor = Internal::ProductOption<MatrixType1, MatrixType2>::Major;
        constexpr static bool isAnyMajor = defaultMajor == MatrixOption::AnyMajor;
        using TargetType = LValueMatrix<OtherDerived>;

        if constexpr (isAnyMajor) {
            for (size_t i = 0; i < target.getMaxMajor(); ++i) {
                for (size_t j = 0; j < target.getMaxMinor(); ++j) {
                    const size_t r = MatrixOption::rowFromMajorMinor<TargetType>(i, j);
                    const size_t c = MatrixOption::columnFromMajorMinor<TargetType>(i, j);
                    target.refFromMajorMinor(i, j) = calc(r, c);
                }
            }
        }
        else {
            for (size_t i = 0; i < (defaultMajor == MatrixOption::Column ? getColumn() : getRow()); ++i) {
                for (size_t j = 0; j < (defaultMajor == MatrixOption::Column ?  getRow() : getColumn()); ++j) {
                    const size_t r = MatrixOption::rowFromMajorMinor<DefaultType>(i, j);
                    const size_t c = MatrixOption::columnFromMajorMinor<DefaultType>(i, j);
                    target(r, c) = calc(r, c);
                }
            }
        }

        constexpr bool isContinuous = std::is_base_of<ContinuousMatrix<OtherDerived>, OtherDerived>::value;
        if constexpr (isContinuous && isReverseDiff)
            target.getDerived().makeContinuous();
    }

    template<class T1, class T2>
    typename MatrixProduct<T1, T2>::ScalarType MatrixProduct<T1, T2>::calc(size_t row, size_t column) const {
        ScalarType result(0);
        for (size_t i = 0; i < mat1.getColumn(); ++i)
            result += ScalarType(mat1.calc(row, i)) * ScalarType(mat2.calc(i, column));
        return result;
    }

    template<class MatrixType1, class MatrixType2>
    [[nodiscard]] inline typename std::enable_if<(MatrixType1::ColumnAtCompile != 1 && MatrixType2::ColumnAtCompile != 1) ||
                                                 (MatrixType1::ColumnAtCompile == 1 && MatrixType2::ColumnAtCompile == 1),
                                                  MatrixProduct<MatrixType1, MatrixType2>>::type
    operator*(const RValueMatrix<MatrixType1>& mat1, const RValueMatrix<MatrixType2>& mat2) noexcept {
        assert(mat1.getColumn() == mat2.getRow());
        return MatrixProduct(mat1, mat2);
    }
}

namespace Physica {
    template<class MatrixType1, class MatrixType2>
    class Traits<MatrixProduct<MatrixType1, MatrixType2>> {
        static_assert(MatrixType1::ColumnAtCompile == MatrixType2::RowAtCompile ||
                      MatrixType1::ColumnAtCompile == Dynamic ||
                      MatrixType2::RowAtCompile == Dynamic,
                      "[Error]: Row and column do not match in matrix-vector product");
    public:
        using ScalarType = typename Core::Internal::BinaryScalarOpReturnType<typename MatrixType1::ScalarType,
                                                                             typename MatrixType2::ScalarType>::Type;
        constexpr static int Option = MatrixOption::AnyMajor | MatrixOption::AnyStorage;
        constexpr static size_t RowAtCompile = MatrixType1::RowAtCompile;
        constexpr static size_t ColumnAtCompile = MatrixType2::ColumnAtCompile;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColumnAtCompile;
    };
}
