/*
 * Copyright 2024 WeiBo He.
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

#include "MatrixImpl/RValueMatrix.h"

namespace Physica::Core {
    template<class ScalarType, size_t Order = Dynamic>
    class UnitMatrix : public RValueMatrix<UnitMatrix<ScalarType, Order>> {
        using This = UnitMatrix<ScalarType, Order>;
        using Base = RValueMatrix<This>;

        size_t order; // TODO: use void if Order != Dynamic
    public:
        UnitMatrix(size_t order_);
        UnitMatrix(const This&) = default;
        UnitMatrix(This&&) noexcept = default;
        ~UnitMatrix() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return ScalarType(row == col ? 1 : 0); }

        [[nodiscard]] const This& transpose() const noexcept { return *this; }
        [[nodiscard]] const This& conjugate() const noexcept { return *this; }
        [[nodiscard]] const This& hermite() const noexcept { return *this; }
        inline void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return order; }
        [[nodiscard]] __host__ __device__ size_t getColumn() const noexcept { return order; }
    };

    template<class ScalarType, size_t Order>
    UnitMatrix<ScalarType, Order>::UnitMatrix(size_t order_) : order(order_) {
        assert((Order == Dynamic || Order == order) && "[Error]: tparam and param is not consistent");
    }

    template<class ScalarType, size_t Order>
    inline void UnitMatrix<ScalarType, Order>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(order, obj.order);
    }

    template<class ScalarType, size_t Order, class VectorType>
    [[nodiscard]] inline const VectorType&
    operator*([[maybe_unused]] const UnitMatrix<ScalarType, Order>& mat, const RValueVector<VectorType>& vec) noexcept {
        assert(mat.getColumn() == vec.getLength());
        return vec.getDerived();
    }
}

namespace Physica {
    using namespace Core;

    template<class T, size_t Order>
    class Traits<UnitMatrix<T, Order>> {
    public:
        using ScalarType = T;
        constexpr static int Option = MatrixOption::AnyMajor | MatrixOption::AnyStorage;
        constexpr static size_t RowAtCompile = Order;
        constexpr static size_t ColumnAtCompile = Order;
        constexpr static size_t MaxRowAtCompile = Order;
        constexpr static size_t MaxColumnAtCompile = Order;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColumnAtCompile;
        constexpr static size_t MaxSizeAtCompile = MaxRowAtCompile * MaxColumnAtCompile;
    };
}
