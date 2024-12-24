/*
 * Copyright 2023-2024 Weibo He.
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
    template<Vector T, int MatrixMajor, size_t Row, size_t Col>
    class ReshapedVector<T, MatrixMajor, Row, Col> : public RValueMatrix<ReshapedVector<T, MatrixMajor, Row, Col>> {
        using This = ReshapedVector<T, MatrixMajor, Row, Col>;
        using Base = RValueMatrix<This>;
    public:
        using ScalarType = T::ScalarType;
    private:
        const T& v;
        size_t r;
        size_t c;
    public:
        ReshapedVector(const T& v_, size_t r_, size_t c_);
        ReshapedVector(const This&) = default;
        ReshapedVector(This&&) noexcept = default;
        ~ReshapedVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t r1, size_t c1) const {
            assert(r1 < getRow() && c1 < getCol());
            if constexpr (MatrixOption::isColMatrix<This>())
                return v.calc(c1 * getRow() + r1);
            else
                return v.calc(r1 * getCol() + c1);
        }
        [[nodiscard]] size_t getRow() const noexcept {
            if constexpr (Row != Dynamic)
                return Row;
            return r;
        }
        [[nodiscard]] size_t getCol() const noexcept {
            if constexpr (Col != Dynamic)
                return Col;
            return c;
        }
        [[nodiscard]] ScalarType sum() const { return v.sum(); }
    };

    template<Vector T, int MatrixMajor, size_t Row, size_t Col>
    ReshapedVector<T, MatrixMajor, Row, Col>::ReshapedVector(const T& v_, size_t r_, size_t c_)
            : v(v_), r(r_), c(c_) {
        assert(r == Row || Row == Dynamic);
        assert(c == Col || Col == Dynamic);
    }

    template<class Derived>
    template<Matrix M>
    auto RValueVector<Derived>::reshape(const M& mat) const {
        return ReshapedVector<Derived, MatrixOption::getMajor<M>(), M::RowAtCompile, M::ColAtCompile>(Base::getDerived(), mat.getRow(), mat.getCol());
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto RValueVector<Derived>::reshape_col(size_t row, size_t col) const {
        return ReshapedVector<Derived, MatrixOption::Col, Row, Col>(Base::getDerived(), row, col);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto RValueVector<Derived>::reshape_row(size_t row, size_t col) const {
        return ReshapedVector<Derived, MatrixOption::Row, Row, Col>(Base::getDerived(), row, col);
    }
}

namespace Physica {
    template<Vector T, int MatrixMajor, size_t Row, size_t Col>
    class Traits<ReshapedVector<T, MatrixMajor, Row, Col>> {
        static_assert(MatrixMajor == MatrixOption::Col || MatrixMajor == MatrixOption::Row, "[Error]: Invalid major");
    public:
        using ScalarType = T::ScalarType;
        constexpr static int Option = MatrixMajor | MatrixOption::AnyStorage;
        constexpr static size_t RowAtCompile = Row;
        constexpr static size_t ColAtCompile = Col;
        constexpr static size_t SizeAtCompile = T::SizeAtCompile;
    };
}
