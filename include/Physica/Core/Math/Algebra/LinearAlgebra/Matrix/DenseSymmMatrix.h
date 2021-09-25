/*
 * Copyright 2020-2021 WeiBo He.
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

#include "LValueMatrix.h"
#include "DenseMatrixImpl/DenseMatrixStorage/HalfDenseMatrixStorage.h"

namespace Physica::Core {
    template<class ScalarType, size_t Order, size_t MaxOrder = Order> class DenseSymmMatrix;

    namespace Internal {
        template<class T> class Traits;

        template<class T, size_t Order, size_t MaxOrder>
        class Traits<DenseSymmMatrix<T, Order, MaxOrder>> {
        public:
            using ScalarType = T;
            constexpr static int MatrixOption = DenseMatrixOption::AnyMajor | DenseMatrixOption::Element;
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
                          , private Internal::HalfDenseMatrixStorage<DenseSymmMatrix<ScalarType, Order, MaxOrder>> {
        using Base = LValueMatrix<DenseSymmMatrix<ScalarType, Order, MaxOrder>>;
        using Storage = Internal::HalfDenseMatrixStorage<DenseSymmMatrix<ScalarType, Order, MaxOrder>>;
    public:
        using Storage::Storage;
        /* Operators */
        using Base::operator=;
        [[nodiscard]] ScalarType& operator()(size_t row, size_t column) { return Storage::operator[](accessingIndex(row, column)); }
        [[nodiscard]] const ScalarType& operator()(size_t row, size_t column) const { return Storage::operator[](accessingIndex(row, column)); }
        /* Getters */
        using Storage::getRow;
        using Storage::getColumn;
    private:
        [[nodiscard]] size_t accessingIndex(size_t r, size_t c) const noexcept;
    };

    template<class ScalarType, size_t Order, size_t MaxOrder>
    size_t DenseSymmMatrix<ScalarType, Order, MaxOrder>::accessingIndex(size_t r, size_t c) const noexcept {
        const size_t order = Storage::getOrder();
        assert(r < order && c < order);
        const bool exchange = c < r;
        const size_t min = exchange ? c : r;
        const size_t max = exchange ? r : c;
        return (order * 2U - min - 1) * min / 2U + max;
    }
}
