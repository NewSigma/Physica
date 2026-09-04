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
    class ReshapedVector : public RValueMatrix<ReshapedVector<V, MatrixMajor, Row, Col>> {
        using This = ReshapedVector<V, MatrixMajor, Row, Col>;
        using Base = RValueMatrix<This>;
        using MaybeRow = std::conditional_t<Row == Dynamic, size_t, Empty>;
        using MaybeCol = std::conditional_t<Col == Dynamic, size_t, Empty>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        decay_rvalue_t<V> v;
        [[no_unique_address]] MaybeRow r;
        [[no_unique_address]] MaybeCol c;
    public:
        ReshapedVector(V v_, size_t r_, size_t c_);
        ReshapedVector(const This&) = default;
        ReshapedVector(This&&) noexcept = default;
        ~ReshapedVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] T calc(size_t r1, size_t c1) const;

        [[nodiscard]] decltype(auto) values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept;
        [[nodiscard]] size_t getCol() const noexcept;
        [[nodiscard]] size_t getOrder() const noexcept;
        [[nodiscard]] T sum() const { return v.sum(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getRowAtCompile() noexcept { return Row; }
        [[nodiscard]] __host__ __device__ consteval static size_t getColAtCompile() noexcept { return Col; }
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept { return MatrixMajor; }
    };

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    ReshapedVector<V, MatrixMajor, Row, Col>::ReshapedVector(V v_, size_t r_, size_t c_)
            : v(std::forward<V>(v_)), r(r_), c(c_) {
        assert(r_ == Row || Row == Dynamic);
        assert(c_ == Col || Col == Dynamic);
        assert(r_ * c_ == v.getLength());
    }

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    auto ReshapedVector<V, MatrixMajor, Row, Col>::calc(size_t r1, size_t c1) const -> T {
        assert(r1 < getRow() && c1 < getCol());
        if constexpr (MatrixMajor::isColMatrix<This>())
            return v.calc(c1 * getRow() + r1);
        else
            return v.calc(r1 * getCol() + c1);
    }

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    decltype(auto) ReshapedVector<V, MatrixMajor, Row, Col>::values(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), V>(self.v).values().template reshape<MatrixMajor, Row, Col>();
    }

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    size_t ReshapedVector<V, MatrixMajor, Row, Col>::getRow() const noexcept {
        if constexpr (Row != Dynamic)
            return Row;
        else
            return r;
    }

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    size_t ReshapedVector<V, MatrixMajor, Row, Col>::getCol() const noexcept {
        if constexpr (Col != Dynamic)
            return Col;
        else
            return c;
    }

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    size_t ReshapedVector<V, MatrixMajor, Row, Col>::getOrder() const noexcept {
        assert(Base::isSquare() && "[Error]: getOrder() assumes square matrix");
        return getRow();
    }
}

namespace Physica {
    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    class Traits<ReshapedVector<V, MatrixMajor, Row, Col>> {
        static_assert(MatrixMajor == MatrixMajor::Col || MatrixMajor == MatrixMajor::Row, "[Error]: Invalid major");
    public:
        using ScalarType = std::remove_cvref_t<V>::ScalarType;
    };
}
