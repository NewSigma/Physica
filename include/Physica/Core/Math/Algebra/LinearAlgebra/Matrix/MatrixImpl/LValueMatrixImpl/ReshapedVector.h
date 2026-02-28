/*
 * Copyright 2025-2026 Weibo He.
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

#include "../LValueMatrix.h"

namespace Physica {
    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    class RValueReshapedVector;

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    class LValueReshapedVector : public LValueMatrix<LValueReshapedVector<V, MatrixMajor, Row, Col>> {
        using This = LValueReshapedVector<V, MatrixMajor, Row, Col>;
        using Base = LValueMatrix<This>;
    public:
        using typename Base::ScalarType;
    private:
        LazyDestroy<V> v;
        size_t r;
        size_t c;
    public:
        LValueReshapedVector(V&& v_, size_t r_, size_t c_);
        LValueReshapedVector(const This&) = default;
        LValueReshapedVector(This&&) noexcept = default;
        ~LValueReshapedVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        using Base::operator=;
        /* Operations */
        void resize(size_t row, size_t col);
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept;
        [[nodiscard]] size_t getCol() const noexcept;
        [[nodiscard]] auto data_ptr(this auto&&, size_t row, size_t col) noexcept;
        [[nodiscard]] ScalarType sum() const { return v.sum(); }
    };

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    LValueReshapedVector<V, MatrixMajor, Row, Col>::LValueReshapedVector(V&& v_, size_t r_, size_t c_)
            : v(std::forward<V>(v_)), r(r_), c(c_) {
        assert(r == Row || Row == Dynamic);
        assert(c == Col || Col == Dynamic);
        assert(r * c == v.getLength());
    }

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    void LValueReshapedVector<V, MatrixMajor, Row, Col>::resize(
            [[maybe_unused]] size_t row, [[maybe_unused]] size_t col) {
        assert(row == getRow());
        assert(col == getCol());
    }

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    size_t LValueReshapedVector<V, MatrixMajor, Row, Col>::getRow() const noexcept {
        if constexpr (Row != Dynamic)
            return Row;
        return r;
    }

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    size_t LValueReshapedVector<V, MatrixMajor, Row, Col>::getCol() const noexcept {
        if constexpr (Col != Dynamic)
            return Col;
        return c;
    }

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    auto LValueReshapedVector<V, MatrixMajor, Row, Col>::data_ptr(this auto&& self, size_t row, size_t col) noexcept {
        assert(row < self.getRow() && col < self.getCol());
        if constexpr (MatrixMajor::isColMatrix<This>())
            return self.v.data_ptr(col * self.getRow() + row);
        else
            return self.v.data_ptr(row * self.getCol() + col);
    }

    template<class Derived>
    auto LValueVector<Derived>::reshape_like(this auto&& self, const Matrix auto& mat) noexcept {
        using Self = decltype(self);
        using M = std::remove_cvref_t<decltype(mat)>;
        constexpr auto Major = MatrixMajor::getMajor<M>();
        static_assert(Major != MatrixMajor::BothMajor, "[Error]: Cannot infer major from this matrix");
        using ResultType = LValueReshapedVector<Self, MatrixMajor::getMajor<M>(), M::RowAtCompile, M::ColAtCompile>;
        return ResultType(std::forward<Self>(self), mat.getRow(), mat.getCol());
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto LValueVector<Derived>::reshape_col(this auto&& self, size_t row, size_t col) noexcept {
        using Self = decltype(self);
        return LValueReshapedVector<Self, MatrixMajor::Col, Row, Col>(std::forward<Self>(self), row, col);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto LValueVector<Derived>::reshape_row(this auto&& self, size_t row, size_t col) noexcept {
        using Self = decltype(self);
        return LValueReshapedVector<Self, MatrixMajor::Row, Row, Col>(std::forward<Self>(self), row, col);
    }
}

namespace Physica {
    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    class Traits<LValueReshapedVector<V, MatrixMajor, Row, Col>>
            : public Traits<RValueReshapedVector<V, MatrixMajor, Row, Col>> {};
}
