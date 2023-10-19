/*
 * Copyright 2020-2022 WeiBo He.
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

#include "DenseMatrix.h"
#include "DenseMatrixImpl/HalfDenseMatrixStorage.h"

namespace Physica::Core {
    template<class ScalarType, size_t Order = Dynamic, size_t MaxOrder = Order> class DenseSymmMatrix;

    namespace Internal {
        template<class T> class Traits;

        template<class T, size_t Order, size_t MaxOrder>
        class Traits<DenseSymmMatrix<T, Order, MaxOrder>> {
        public:
            using ScalarType = T;
            constexpr static int Option = MatrixOption::AnyMajor | MatrixOption::Element;
            constexpr static size_t RowAtCompile = Order;
            constexpr static size_t ColumnAtCompile = Order;
            constexpr static size_t MaxRowAtCompile = MaxOrder;
            constexpr static size_t MaxColumnAtCompile = MaxOrder;
            constexpr static size_t SizeAtCompile = RowAtCompile * ColumnAtCompile;
            constexpr static size_t MaxSizeAtCompile = MaxRowAtCompile * MaxColumnAtCompile;
        };
    }

    template<class ScalarType, size_t Order, size_t MaxOrder>
    class DenseSymmMatrix : public LValueMatrix<DenseSymmMatrix<ScalarType, Order, MaxOrder>>
                          , private Internal::HalfDenseMatrixStorage<ScalarType, Order, MaxOrder> {
        using Base = LValueMatrix<DenseSymmMatrix<ScalarType, Order, MaxOrder>>;
        using Storage = Internal::HalfDenseMatrixStorage<ScalarType, Order, MaxOrder>;
    public:
        using ColMatrix = DenseSymmMatrix<ScalarType, Order, MaxOrder>;
        using RowMatrix = DenseSymmMatrix<ScalarType, Order, MaxOrder>;
        using RealMatrix = DenseSymmMatrix<typename ScalarType::RealType, Order, MaxOrder>;
    public:
        template<class OtherMatrix>
        DenseSymmMatrix(const RValueMatrix<OtherMatrix>& mat);
        using Storage::Storage;
        DenseSymmMatrix(const DenseSymmMatrix&) = default;
        DenseSymmMatrix(DenseSymmMatrix&&) noexcept = default;
        ~DenseSymmMatrix() = default;
        /* Operators */
        DenseSymmMatrix& operator=(DenseSymmMatrix m) noexcept;
        using Base::operator=;
        [[nodiscard]] ScalarType& operator()(size_t row, size_t column);
        [[nodiscard]] const ScalarType& operator()(size_t row, size_t column) const;
        /* Operations */
        using Storage::resize;
        /* Getters */
        using Storage::getOrder;
        using Storage::getRow;
        using Storage::getColumn;
        /* Helpers */
        void swap(DenseSymmMatrix& m) noexcept;
        /* Static members */
        [[nodiscard]] static DenseSymmMatrix unitMatrix(size_t order);
    };
    /**
     * Assuming mat is a symmetric matrix, if it is not the case, only half of the elements is saved correctly
     */
    template<class ScalarType, size_t Order, size_t MaxOrder>
    template<class OtherMatrix>
    DenseSymmMatrix<ScalarType, Order, MaxOrder>::DenseSymmMatrix(const RValueMatrix<OtherMatrix>& mat)
            : DenseSymmMatrix(mat.getRow()) {
        assert(mat.getRow() == mat.getColumn());
        for (size_t i = 0; i < mat.getRow(); ++i)
            for (size_t j = i; j < mat.getRow(); ++j)
                operator()(i, j) = mat.calc(i, j);
    }

    template<class ScalarType, size_t Order, size_t MaxOrder>
    DenseSymmMatrix<ScalarType, Order, MaxOrder>&
    DenseSymmMatrix<ScalarType, Order, MaxOrder>::operator=(DenseSymmMatrix m) noexcept {
        swap(m);
        return *this;
    }

    template<class ScalarType, size_t Order, size_t MaxOrder>
    ScalarType& DenseSymmMatrix<ScalarType, Order, MaxOrder>::operator()(size_t row, size_t column) {
        const size_t index = Storage::accessingIndex(row, column);
        return Storage::operator[](index);
    }

    template<class ScalarType, size_t Order, size_t MaxOrder>
    const ScalarType& DenseSymmMatrix<ScalarType, Order, MaxOrder>::operator()(size_t row, size_t column) const {
        const size_t index = Storage::accessingIndex(row, column);
        return Storage::operator[](index);
    }

    template<class ScalarType, size_t Order, size_t MaxOrder>
    void DenseSymmMatrix<ScalarType, Order, MaxOrder>::swap(DenseSymmMatrix& m) noexcept {
        assert(this != &m && "[Error]: Self swap is likely a bug");
        Storage::swap(m);
    }

    template<class ScalarType, size_t Order, size_t MaxOrder>
    DenseSymmMatrix<ScalarType, Order, MaxOrder> DenseSymmMatrix<ScalarType, Order, MaxOrder>::unitMatrix(size_t order) {
        DenseSymmMatrix<ScalarType, Order, MaxOrder> result(order);
        result.toUnitMatrix();
        return result;
    }

    template<class ScalarType, size_t Order, size_t MaxOrder>
    inline void swap(Physica::Core::DenseSymmMatrix<ScalarType, Order, MaxOrder>& m1,
                     Physica::Core::DenseSymmMatrix<ScalarType, Order, MaxOrder>& m2) noexcept {
        m1.swap(m2);
    }
}
