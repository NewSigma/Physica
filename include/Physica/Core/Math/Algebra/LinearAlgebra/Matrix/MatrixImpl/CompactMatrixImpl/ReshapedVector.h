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

#include "../CompactMatrix.h"

namespace Physica {
    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    class CompactReshapedVector : public CompactMatrix<CompactReshapedVector<V, MatrixMajor, Row, Col>> {
        using This = CompactReshapedVector<V, MatrixMajor, Row, Col>;
        using Base = LValueMatrix<This>;
        using MaybeRow = std::conditional_t<Row == Dynamic, size_t, Empty>;
        using MaybeCol = std::conditional_t<Col == Dynamic, size_t, Empty>;
    protected:
        using typename Base::T;
    private:
        decay_rvalue_t<V> v;
        [[no_unique_address]] MaybeRow r;
        [[no_unique_address]] MaybeCol c;
    public:
        CompactReshapedVector(V&& v_, size_t r_, size_t c_);
        CompactReshapedVector(const This&) = default;
        CompactReshapedVector(This&&) noexcept = default;
        ~CompactReshapedVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        using Base::operator=;
        /* Operations */
        void resize(size_t row, size_t col);

        [[nodiscard]] T sum() const { return v.sum(); }
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept;
        [[nodiscard]] size_t getCol() const noexcept;
        [[nodiscard]] size_t getOrder() const noexcept;
        [[nodiscard]] auto data(this auto&&) noexcept;
        [[nodiscard]] auto data_ptr(this auto&&, size_t row, size_t col) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept { return MatrixMajor; }
        /* Friends */
        friend class device_obj<This>;
    };

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    CompactReshapedVector<V, MatrixMajor, Row, Col>::CompactReshapedVector(V&& v_, size_t r_, size_t c_)
            : v(std::forward<V>(v_)), r(r_), c(c_) {
        assert(r == Row || Row == Dynamic);
        assert(c == Col || Col == Dynamic);
        assert(r * c == v.getLength());
    }

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    void CompactReshapedVector<V, MatrixMajor, Row, Col>::resize(
            [[maybe_unused]] size_t row, [[maybe_unused]] size_t col) {
        assert(row == getRow());
        assert(col == getCol());
    }

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    size_t CompactReshapedVector<V, MatrixMajor, Row, Col>::getRow() const noexcept {
        if constexpr (Row != Dynamic)
            return Row;
        else
            return r;
    }

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    size_t CompactReshapedVector<V, MatrixMajor, Row, Col>::getCol() const noexcept {
        if constexpr (Col != Dynamic)
            return Col;
        else
            return c;
    }

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    auto CompactReshapedVector<V, MatrixMajor, Row, Col>::data(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self).v.data();
    }

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    auto CompactReshapedVector<V, MatrixMajor, Row, Col>::data_ptr(this auto&& self, size_t row, size_t col) noexcept {
        if constexpr (self.isColMatrix())
            return self.data() + col * self.getRow() + row;
        else
            return self.data() + row * self.getCol() + col;
    }

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    size_t CompactReshapedVector<V, MatrixMajor, Row, Col>::getOrder() const noexcept {
        assert(Base::isSquare() && "[Error]: getOrder() assumes square matrix");
        return getRow();
    }
}

namespace Physica {
    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    class Traits<CompactReshapedVector<V, MatrixMajor, Row, Col>> : public Traits<RValueReshapedVector<V, MatrixMajor, Row, Col>> {};
}
