/*
 * Copyright 2020-2024 WeiBo He.
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

#include <memory>
#include "DenseMatrixImpl/DenseMatrixDim.h"
#include "DenseMatrixImpl/DenseMatrixStorage.h"
#include "InverseMatrix.h"
#include "MatrixDecomposition/LUDecomposition.h"
#include "MatrixImpl/ContinuousMatrix.h"
#include "MatrixProduct/MatrixProduct.h"

namespace Physica::Core {
    /**
     * \class DenseMatrix
     * A matrix can be either fixed matrix, which have its max size defined,
     * or dynamic matrix, whose size is dynamically changed.
     *
     * \tparam Option
     * Option is combinations of \enum MatrixOption
     */
    template<class T,
            int Option = MatrixOption::Column | MatrixOption::Vector,
            size_t Row = Dynamic,
            size_t Column = Dynamic,
            size_t MaxRow = Row,
            size_t MaxColumn = Column,
            class Allocator = Utils::HostAllocator<T>>
    class DenseMatrix : public ContinuousMatrix<DenseMatrix<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>>
                      , public DenseMatrixStorage<DenseMatrix<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>, Option>
                      , public DenseMatrixDim<DenseMatrix<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>, Row, Column, MaxRow, MaxColumn> {
        using This = DenseMatrix<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>;
        using Base = ContinuousMatrix<This>;
        using Storage = DenseMatrixStorage<This, Option>;
        using Dim = DenseMatrixDim<This, Row, Column, MaxRow, MaxColumn>;
        using InitializerType = typename Storage::InitializerType;
        using Base::isReverseDiff;
    public:
        using typename Base::PlainScalar;
        using typename Base::ScalarType;
        using device_obj_type = device_obj<This>;
        using ColMatrix = DenseMatrix<T, MatrixOption::getStorage<DenseMatrix>() | MatrixOption::Column, Row, Column, MaxRow, MaxColumn>;
        using RowMatrix = DenseMatrix<T, MatrixOption::getStorage<DenseMatrix>() | MatrixOption::Row, Row, Column, MaxRow, MaxColumn>;
        using RealMatrix = DenseMatrix<typename T::RealType, Option, Row, Column, MaxRow, MaxColumn>;
    public:
        DenseMatrix() = default;
        DenseMatrix(size_t row, size_t column);
        DenseMatrix(size_t row, size_t column, T value);
        DenseMatrix(std::initializer_list<InitializerType> list);
        template<class OtherMatrix>
        DenseMatrix(const RValueMatrix<OtherMatrix>& mat);
        template<class VectorType>
        DenseMatrix(const RValueVector<VectorType>& mat);
        template<class MatrixIn>
        DenseMatrix(LUDecomposition<MatrixIn> lu);
        DenseMatrix(const This& m);
        DenseMatrix(This&& m) noexcept;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        using Base::operator=;
        using Storage::operator();
        /* Operations */
        size_t completePivoting(size_t column);
        size_t partialPivoting(size_t column);

        inline void resize(size_t row, size_t column);
        [[nodiscard]] DenseMatrix copy() const;
        [[nodiscard]] inline device_obj<This> toDevice() const;
        void toDevice(device_obj<This>& obj) const;
        void swap(DenseMatrix& __restrict m) noexcept;

        using Base::random_any;
        using Base::random_normal;
        using Base::random_uniform;
        /* Getters */
        using Dim::getColumn;
        using Dim::getRow;
        using Storage::data_ptr;
        /* Static members */
        [[nodiscard]] static DenseMatrix Zeros(size_t rank) { return DenseMatrix(rank, rank, T(0)); }
        [[nodiscard]] static DenseMatrix Zeros(size_t row, size_t column) { return DenseMatrix(row, column, T(0)); }
        [[nodiscard]] static DenseMatrix unitMatrix(size_t order);
        template<class RandomGenerator>
        [[nodiscard]] static DenseMatrix random_uniform(size_t order, RandomGenerator& gen) { return random_uniform(order, order, gen); }
        template<class RandomGenerator>
        [[nodiscard]] inline static DenseMatrix random_uniform(size_t row, size_t column, RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] inline static DenseMatrix random_normal(size_t row, size_t column, RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        [[nodiscard]] inline static DenseMatrix random_any(size_t row, size_t column, Distribution& dist, RandomGenerator& gen);
        template<class VectorType>
        [[nodiscard]] static std::pair<DenseMatrix, DenseMatrix> meshgrid(const LValueVector<VectorType>& vecInCols, const LValueVector<VectorType>& vecInRows);
    private:
        DenseMatrix(Storage storage, size_t row, size_t column) : Storage(std::move(storage)), Dim(row, column) {}
        friend class device_obj<This>;
    };

    template<class T, int Option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    std::istream& operator>>(std::istream& is, DenseMatrix<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>& mat);
}

namespace Physica {
    template<class T, int Op, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    class Traits<Core::DenseMatrix<T, Op, Row, Column, MaxRow, MaxColumn, Allocator>> {
        static_assert(MaxRow * MaxColumn * sizeof(T) <= 2048, "[Warn]: It is suggested declare large fixed size matrix as dynamic matrix");
    public:
        using ScalarType = T;
        constexpr static int Option = Op;
        constexpr static size_t RowAtCompile = Row;
        constexpr static size_t ColumnAtCompile = Column;
        constexpr static size_t MaxRowAtCompile = MaxRow;
        constexpr static size_t MaxColumnAtCompile = MaxColumn;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColumnAtCompile;
        constexpr static size_t MaxSizeAtCompile = MaxRowAtCompile * MaxColumnAtCompile;
        using AllocatorType = Allocator;
    };
}

namespace std {
    template<class T, int Option, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator>
    inline void swap(
            Physica::Core::DenseMatrix<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>& __restrict m1,
            Physica::Core::DenseMatrix<T, Option, Row, Column, MaxRow, MaxColumn, Allocator>& __restrict m2) noexcept {
        m1.swap(m2);
    }
}

#include "DenseMatrixImpl/DenseMatrixImpl.h"
