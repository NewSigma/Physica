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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/CompactVector.h"
#include "Physica/Core/Utils/Empty.h"

namespace Physica {
    template<class MatrixType, size_t Row = Dynamic, size_t Col = Dynamic> class CompactMatrixBlock;

    template<Matrix M, size_t Row, size_t Col>
    class CompactMatrixBlock<M, Row, Col> : public CompactVector<CompactMatrixBlock<M, Row, Col>> {
        using This = CompactMatrixBlock<M, Row, Col>;
        using Base = CompactVector<This>;
        using MaybeRowCount = std::conditional_t<Row == Dynamic, size_t, Empty>;
        using MaybeColCount = std::conditional_t<Col == Dynamic, size_t, Empty>;
    private:
        decay_rvalue_t<M> mat;
        size_t fromRow;
        [[no_unique_address]] MaybeRowCount rowCount;
        size_t fromCol;
        [[no_unique_address]] MaybeColCount colCount;
    public:
        constexpr CompactMatrixBlock(M&& mat_, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount);
        CompactMatrixBlock(const This&) = default;
        CompactMatrixBlock(This&&) noexcept = default;
        ~CompactMatrixBlock() = default;
        /* Operators */
        This& operator=(const This& m);
        This& operator=(This&& m) noexcept;
        using Base::operator=;
        /* Operations */
        using Base::resize;
        constexpr void resize(size_t length);

        [[nodiscard]] auto&& row(this auto&&, size_t r) noexcept requires(Row == 1);
        [[nodiscard]] auto&& col(this auto&&, size_t c) noexcept requires(Col == 1);

        [[nodiscard]] auto values(this auto&&) noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] auto grads(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] constexpr size_t getLength() const noexcept;
        [[nodiscard]] constexpr auto data_handle(this auto&& self) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept;
    };

    template<Matrix M, size_t Row, size_t Col>
    constexpr CompactMatrixBlock<M, Row, Col>::CompactMatrixBlock(M&& mat_, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount)
            : mat(std::forward<M>(mat_))
            , fromRow(fromRow)
            , rowCount(rowCount)
            , fromCol(fromCol)
            , colCount(colCount) {
        assert(fromRow < mat.getRow());
        assert(fromCol < mat.getCol());
        assert(fromRow + rowCount <= mat.getRow());
        assert(fromCol + colCount <= mat.getCol());
        if constexpr (Row != Dynamic)
            assert(Row == rowCount);
        if constexpr (Col != Dynamic)
            assert(Col == colCount);
    }

    template<Matrix M, size_t Row, size_t Col>
    auto CompactMatrixBlock<M, Row, Col>::operator=(const This& m) -> This& {
        Base::operator=(m);
        return *this;
    }

    template<Matrix M, size_t Row, size_t Col>
    auto CompactMatrixBlock<M, Row, Col>::operator=(This&& m) noexcept -> This& {
        return *this = m;
    }

    template<Matrix M, size_t Row, size_t Col>
    constexpr void CompactMatrixBlock<M, Row, Col>::resize([[maybe_unused]] size_t length) {
        assert(length == getLength());
    }

    template<Matrix M, size_t Row, size_t Col>
    auto&& CompactMatrixBlock<M, Row, Col>::row(this auto&& self, [[maybe_unused]] size_t r) noexcept requires(Row == 1) {
        assert(r == 0);
        return std::forward<decltype(self)>(self);
    }

    template<Matrix M, size_t Row, size_t Col>
    auto&& CompactMatrixBlock<M, Row, Col>::col(this auto&& self, [[maybe_unused]] size_t c) noexcept requires(Col == 1) {
        assert(c == 0);
        return std::forward<decltype(self)>(self);
    }

    template<Matrix M, size_t Row, size_t Col>
    auto CompactMatrixBlock<M, Row, Col>::values(this auto&& self) noexcept {
        auto&& v = propagate_rvalue_reference<decltype(self), M>(self.mat).values();
        using M1 = decltype(v);
        return CompactMatrixBlock<M1, Row, Col>(std::forward<M1>(v), self.fromRow, Row == 1 ? 1 : self.getLength(), self.fromCol, Col == 1 ? 1 : self.getLength());
    }

    template<Matrix M, size_t Row, size_t Col>
    template<int GradOrder>
    auto CompactMatrixBlock<M, Row, Col>::grads(this auto&& self) noexcept {
        auto&& g = propagate_rvalue_reference<decltype(self), M>(self.mat).template grads<GradOrder>();
        using M1 = decltype(g);
        return CompactMatrixBlock<M1, Row, Col>(std::forward<M1>(g), self.fromRow, Row == 1 ? 1 : self.getLength(), self.fromCol, Col == 1 ? 1 : self.getLength());
    }

    template<Matrix M, size_t Row, size_t Col>
    constexpr size_t CompactMatrixBlock<M, Row, Col>::getLength() const noexcept {
        if constexpr (Row == 1) {
            if constexpr (Col == Dynamic)
                return colCount;
            else
                return Col;
        }
        else {
            if constexpr (Row == Dynamic)
                return rowCount;
            else
                return Row;
        }
    }

    template<Matrix M, size_t Row, size_t Col>
    constexpr auto CompactMatrixBlock<M, Row, Col>::data_handle(this auto&& self) noexcept {
        return self.mat.data_ptr(self.fromRow, self.fromCol);
    }

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ consteval size_t CompactMatrixBlock<M, Row, Col>::getSizeAtCompile() noexcept {
        return Row == 1 ? Col : Row;
    }
}

namespace Physica {
    template<Matrix M, size_t Row, size_t Col>
    class Traits<CompactMatrixBlock<M, Row, Col>> {
    public:
        using ScalarType = std::remove_cvref_t<M>::ScalarType;
    };
}
