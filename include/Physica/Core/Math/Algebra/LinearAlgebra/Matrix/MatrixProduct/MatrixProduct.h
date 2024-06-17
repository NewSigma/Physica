/*
 * Copyright 2021-2023 WeiBo He.
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
    template<class Derived> class ContinuousMatrix;
    template<class MatrixType1, class MatrixType2> class MatrixProduct;
    template<class VectorType, class MatrixType> class VectorMatrixProduct;
    template<class MatrixType, class VectorType> class MatrixVectorProduct;

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
            constexpr static int Option = MatrixOption::AnyMajor | MatrixOption::AnyStorage;
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
                      "[Error]: Row and column do not match in matrix-vector product");
        using This = MatrixProduct<MatrixType1, MatrixType2>;
    public:
        using Base = RValueMatrix<This>;
        using Base::isReverseDiff;
        using typename Base::ScalarType;
        using DefaultType = DenseMatrix<ScalarType,
                                        Internal::ProductOption<MatrixType1, MatrixType2>::Option,
                                        Base::RowAtCompile,
                                        Base::ColumnAtCompile,
                                        Base::MaxRowAtCompile,
                                        Base::MaxColumnAtCompile,
                                        Utils::HostAllocator<ScalarType>>;
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

    template<class VectorType, class MatrixType>
    class VectorMatrixProduct : public RValueMatrix<VectorMatrixProduct<VectorType, MatrixType>> {
        static_assert(MatrixType::RowAtCompile == 1 || MatrixType::RowAtCompile == Dynamic,
                      "Row and column do not match in matrix product");
        using This = VectorMatrixProduct<VectorType, MatrixType>;
    public:
        using Base = RValueMatrix<This>;
        using Base::isReverseDiff;
        using typename Base::ScalarType;
    private:
        const VectorType& vec;
        const MatrixType& mat;
    public:
        VectorMatrixProduct(const RValueVector<VectorType>& vec_, const RValueMatrix<MatrixType>& mat_)
                : vec(vec_.getDerived()), mat(mat_.getDerived()) {
            assert(mat.getRow() == 1);
        }
        VectorMatrixProduct(const This&) = delete;
        VectorMatrixProduct(This&&) noexcept = delete;
        ~VectorMatrixProduct() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t column) const;
        [[nodiscard]] __host__ __device__ size_t getRow() const { return vec.getLength(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return mat.getColumn(); }
        [[nodiscard]] const VectorType& getLHS() const noexcept { return vec; }
        [[nodiscard]] const MatrixType& getRHS() const noexcept { return mat; }
    };

    template<class MatrixType, class VectorType>
    class MatrixVectorProduct : public RValueVector<MatrixVectorProduct<MatrixType, VectorType>> {
        static_assert(MatrixType::ColumnAtCompile == VectorType::SizeAtCompile ||
                      MatrixType::ColumnAtCompile == Dynamic ||
                      VectorType::SizeAtCompile == Dynamic,
                      "Row and column do not match in matrix product");
        using This = MatrixVectorProduct<MatrixType, VectorType>;
    public:
        using Base = RValueVector<This>;
        using typename Base::ScalarType;
    private:
        const MatrixType& mat;
        const VectorType& vec;
    public:
        MatrixVectorProduct(const RValueMatrix<MatrixType>& mat_, const RValueVector<VectorType>& vec_)
                : mat(mat_.getDerived()), vec(vec_.getDerived()) {
            assert(mat.getColumn() == vec.getLength());
        }
        MatrixVectorProduct(const This&) = delete;
        MatrixVectorProduct(This&&) noexcept = delete;
        ~MatrixVectorProduct() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<class OtherDerived>
        inline void assignTo(LValueVector<OtherDerived>& target) const;
        /* Getters */
        [[nodiscard]] inline ScalarType calc(size_t index) const;
        [[nodiscard]] __host__ __device__ size_t getLength() const { return mat.getRow(); }
        [[nodiscard]] const MatrixType& getLHS() const noexcept { return mat; }
        [[nodiscard]] const VectorType& getRHS() const noexcept { return vec; }
    };

    template<class MatrixType1, class MatrixType2>
    [[nodiscard]] inline typename std::enable_if<(MatrixType1::ColumnAtCompile != 1 && MatrixType2::ColumnAtCompile != 1) ||
                                                 (MatrixType1::ColumnAtCompile == 1 && MatrixType2::ColumnAtCompile == 1),
                                                  MatrixProduct<MatrixType1, MatrixType2>>::type
    operator*(const RValueMatrix<MatrixType1>& mat1, const RValueMatrix<MatrixType2>& mat2) noexcept {
        assert(mat1.getColumn() == mat2.getRow());
        return MatrixProduct(mat1, mat2);
    }
    /**
     * \note Here we force the row of \param mat is 1, because in Physica vectors are naturally column vectors.
     * To compute row vector * matrix, users should converted it to matrix^T * column vector.
     */
    template<class VectorType, class MatrixType>
    [[nodiscard]] inline typename std::enable_if<MatrixType::RowAtCompile == 1, VectorMatrixProduct<VectorType, MatrixType>>::type
    operator*(const RValueVector<VectorType>& vec, const RValueMatrix<MatrixType>& mat) noexcept {
        assert(mat.getRow() == 1);
        return VectorMatrixProduct(vec, mat);
    }

    template<class MatrixType, class VectorType>
    [[nodiscard]] inline typename std::enable_if<MatrixType::RowAtCompile != 1, MatrixVectorProduct<MatrixType, VectorType>>::type
    operator*(const RValueMatrix<MatrixType>& mat, const RValueVector<VectorType>& vec) noexcept {
        assert(mat.getColumn() == vec.getLength());
        return MatrixVectorProduct(mat.getDerived(), vec.getDerived());
    }

    template<class MatrixType, class VectorType>
    [[nodiscard]] inline typename std::enable_if<MatrixType::RowAtCompile == 1 && MatrixType::ColumnAtCompile == 1,
                                                 typename Internal::BinaryScalarOpReturnType<typename MatrixType::ScalarType,
                                                                                             typename VectorType::ScalarType>::Type>::type
    operator*(const RValueMatrix<MatrixType>& mat, const RValueVector<VectorType>& vec) {
        assert(mat.getColumn() == vec.getLength());
        return mat.row(0) * vec;
    }
}

#include "MatrixProductImpl.h"
