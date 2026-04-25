/*
 * Copyright 2023-2026 Weibo He.
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

#include "../RValueMatrix.h"

namespace Physica {
    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    class RValueReshapedVector : public RValueMatrix<RValueReshapedVector<V, MatrixMajor, Row, Col>> {
        using This = RValueReshapedVector<V, MatrixMajor, Row, Col>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        LazyDestroy<V> v;
        size_t r;
        size_t c;
    public:
        RValueReshapedVector(V v_, size_t r_, size_t c_);
        RValueReshapedVector(const This&) = default;
        RValueReshapedVector(This&&) noexcept = default;
        ~RValueReshapedVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] T calc(size_t r1, size_t c1) const;

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept;
        [[nodiscard]] size_t getCol() const noexcept;
        [[nodiscard]] T sum() const { return v.sum(); }
    };

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    RValueReshapedVector<V, MatrixMajor, Row, Col>::RValueReshapedVector(V v_, size_t r_, size_t c_)
            : v(std::forward<V>(v_)), r(r_), c(c_) {
        assert(r == Row || Row == Dynamic);
        assert(c == Col || Col == Dynamic);
        assert(r * c == v.getLength());
    }

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    auto RValueReshapedVector<V, MatrixMajor, Row, Col>::calc(size_t r1, size_t c1) const -> T {
        assert(r1 < getRow() && c1 < getCol());
        if constexpr (MatrixMajor::isColMatrix<This>())
            return v.calc(c1 * getRow() + r1);
        else
            return v.calc(r1 * getCol() + c1);
    }

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    auto RValueReshapedVector<V, MatrixMajor, Row, Col>::values(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), V>(self.v).values().template reshape<MatrixMajor, Row, Col>();
    }

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    size_t RValueReshapedVector<V, MatrixMajor, Row, Col>::getRow() const noexcept {
        if constexpr (Row != Dynamic)
            return Row;
        return r;
    }

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    size_t RValueReshapedVector<V, MatrixMajor, Row, Col>::getCol() const noexcept {
        if constexpr (Col != Dynamic)
            return Col;
        return c;
    }
}

namespace Physica {
    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    class Traits<RValueReshapedVector<V, MatrixMajor, Row, Col>> {
        static_assert(MatrixMajor == MatrixMajor::Col || MatrixMajor == MatrixMajor::Row, "[Error]: Invalid major");
    public:
        using ScalarType = std::remove_cvref_t<V>::ScalarType;
        constexpr static int Major = MatrixMajor;
        constexpr static size_t RowAtCompile = Row;
        constexpr static size_t ColAtCompile = Col;
        constexpr static size_t SizeAtCompile = std::remove_cvref_t<V>::getSizeAtCompile();
    };
}
