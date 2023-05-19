/*
 * Copyright 2023 WeiBo He.
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

#include "RValueMatrix.h"

namespace Physica::Core {
    template<class VectorType, int MatrixMajor, size_t Row, size_t Column> class ReshapedVector;

    namespace Internal {
        template<class T> class Traits;

        template<class VectorType, int MatrixMajor, size_t Row, size_t Column>
        class Traits<ReshapedVector<VectorType, MatrixMajor, Row, Column>> {
            static_assert(MatrixMajor == MatrixOption::Column || MatrixMajor == MatrixOption::Row, "[Error]: Invalid major");
        public:
            using ScalarType = typename VectorType::ScalarType;
            constexpr static int Option = MatrixMajor | MatrixOption::AnyStorage;
            constexpr static size_t RowAtCompile = Row;
            constexpr static size_t ColumnAtCompile = Column;
            constexpr static size_t MaxRowAtCompile = RowAtCompile;
            constexpr static size_t MaxColumnAtCompile = ColumnAtCompile;
            constexpr static size_t SizeAtCompile = VectorType::SizeAtCompile;
            constexpr static size_t MaxSizeAtCompile = SizeAtCompile;
        };
    }

    template<class VectorType, int MatrixMajor, size_t Row, size_t Column>
    class ReshapedVector : public RValueMatrix<ReshapedVector<VectorType, MatrixMajor, Row, Column>> {
        using This = ReshapedVector<VectorType, MatrixMajor, Row, Column>;
        using Base = RValueMatrix<This>;
        using ScalarType = typename VectorType::ScalarType;
    public:
        const VectorType& v;
        size_t row;
        size_t column;
    public:
        ReshapedVector(const VectorType& v_, size_t row, size_t column);
        ReshapedVector(const ReshapedVector&) = default;
        ReshapedVector(ReshapedVector&&) noexcept = default;
        ~ReshapedVector() = default;
        /* Operators */
        ReshapedVector& operator=(const ReshapedVector&) = delete;
        ReshapedVector& operator=(ReshapedVector&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t r, size_t c) const {
            assert(r < getRow() && c < getColumn());
            if constexpr (MatrixOption::isColumnMatrix<This>())
                return v.calc(c * getRow() + r);
            else
                return v.calc(r * getColumn() + c);
        }
        [[nodiscard]] size_t getRow() const noexcept {
            if constexpr (Row != Dynamic)
                return Row;
            return row;
        }
        [[nodiscard]] size_t getColumn() const noexcept {
            if constexpr (Column != Dynamic)
                return Column;
            return column;
        }
        [[nodiscard]] ScalarType sum() const { return v.sum(); }
    };

    template<class VectorType, int MatrixMajor, size_t Row, size_t Column>
    ReshapedVector<VectorType, MatrixMajor, Row, Column>::ReshapedVector(const VectorType& v_, size_t row_, size_t column_)
            : v(v_), row(row_), column(column_) {}

    template<class Derived>
    template<class OtherDerived>
    ReshapedVector<Derived, MatrixOption::getMajor<OtherDerived>(), OtherDerived::RowAtCompile, OtherDerived::ColumnAtCompile>
    RValueVector<Derived>::reshape(const RValueMatrix<OtherDerived>& mat) const {
        return {Base::getDerived(), mat.getRow(), mat.getColumn()};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    ReshapedVector<Derived, Row, Column, MatrixOption::Column> RValueVector<Derived>::reshape_col(size_t row, size_t col) const {
        return {Base::getDerived(), row, col};
    }

    template<class Derived>
    template<size_t Row, size_t Column>
    ReshapedVector<Derived, Row, Column, MatrixOption::Row> RValueVector<Derived>::reshape_row(size_t row, size_t col) const {
        return {Base::getDerived(), row, col};
    }
}
