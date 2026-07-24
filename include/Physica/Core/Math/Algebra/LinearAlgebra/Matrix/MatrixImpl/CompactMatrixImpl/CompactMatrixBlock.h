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

#include "../CompactMatrix.h"
#include "Physica/Core/Utils/Empty.h"

namespace Physica {
    template<class MatrixType, size_t Row = Dynamic, size_t Col = Dynamic> class CompactMatrixBlock;

    template<Matrix M, size_t Row, size_t Col> requires(Row == 1 || Col == 1)
    class CompactMatrixBlock<M, Row, Col> : public CompactVector<CompactMatrixBlock<M, Row, Col>> {
        static_assert((Row == 1 && MatrixMajor::isRowMatrix<M>()) || (Col == 1 && MatrixMajor::isColMatrix<M>()) || (Row == 1 || Col == 1),
                      "[Error]: The matrix block is not compact");
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
        CompactMatrixBlock(M&& mat_, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount);
        CompactMatrixBlock(const This&) = default;
        CompactMatrixBlock(This&&) noexcept = default;
        ~CompactMatrixBlock() = default;
        /* Operators */
        This& operator=(const This& m);
        This& operator=(This&& m) noexcept;
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize(size_t length);

        [[nodiscard]] auto&& row(this auto&&, size_t r) noexcept requires(Row == 1);
        [[nodiscard]] auto&& col(this auto&&, size_t c) noexcept requires(Col == 1);

        [[nodiscard]] auto values(this auto&&) noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] auto grads(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] constexpr size_t getLength() const noexcept;
        [[nodiscard]] auto data(this auto&& self) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept;
    };

    template<Matrix M, size_t Row, size_t Col> requires(Row == 1 || Col == 1)
    CompactMatrixBlock<M, Row, Col>::CompactMatrixBlock(M&& mat_, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount)
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

    template<Matrix M, size_t Row, size_t Col> requires(Row == 1 || Col == 1)
    auto CompactMatrixBlock<M, Row, Col>::operator=(const This& m) -> This& {
        Base::operator=(m);
        return *this;
    }

    template<Matrix M, size_t Row, size_t Col> requires(Row == 1 || Col == 1)
    auto CompactMatrixBlock<M, Row, Col>::operator=(This&& m) noexcept -> This& {
        return *this = m;
    }

    template<Matrix M, size_t Row, size_t Col> requires(Row == 1 || Col == 1)
    void CompactMatrixBlock<M, Row, Col>::resize([[maybe_unused]] size_t length) {
        assert(length == getLength());
    }

    template<Matrix M, size_t Row, size_t Col> requires(Row == 1 || Col == 1)
    auto&& CompactMatrixBlock<M, Row, Col>::row(this auto&& self, [[maybe_unused]] size_t r) noexcept requires(Row == 1) {
        assert(r == 0);
        return std::forward<decltype(self)>(self);
    }

    template<Matrix M, size_t Row, size_t Col> requires(Row == 1 || Col == 1)
    auto&& CompactMatrixBlock<M, Row, Col>::col(this auto&& self, [[maybe_unused]] size_t c) noexcept requires(Col == 1) {
        assert(c == 0);
        return std::forward<decltype(self)>(self);
    }

    template<Matrix M, size_t Row, size_t Col> requires(Row == 1 || Col == 1)
    auto CompactMatrixBlock<M, Row, Col>::values(this auto&& self) noexcept {
        auto&& v = propagate_rvalue_reference<decltype(self), M>(self.mat).values();
        using M1 = decltype(v);
        return CompactMatrixBlock<M1, Row, Col>(std::forward<M1>(v), self.fromRow, Row == 1 ? 1 : self.getLength(), self.fromCol, Col == 1 ? 1 : self.getLength());
    }

    template<Matrix M, size_t Row, size_t Col> requires(Row == 1 || Col == 1)
    template<int GradOrder>
    auto CompactMatrixBlock<M, Row, Col>::grads(this auto&& self) noexcept {
        auto&& g = propagate_rvalue_reference<decltype(self), M>(self.mat).template grads<GradOrder>();
        using M1 = decltype(g);
        return CompactMatrixBlock<M1, Row, Col>(std::forward<M1>(g), self.fromRow, Row == 1 ? 1 : self.getLength(), self.fromCol, Col == 1 ? 1 : self.getLength());
    }

    template<Matrix M, size_t Row, size_t Col> requires(Row == 1 || Col == 1)
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

    template<Matrix M, size_t Row, size_t Col> requires(Row == 1 || Col == 1)
    auto CompactMatrixBlock<M, Row, Col>::data(this auto&& self) noexcept {
        return self.mat.data_ptr(self.fromRow, self.fromCol);
    }

    template<Matrix M, size_t Row, size_t Col> requires(Row == 1 || Col == 1)
    __host__ __device__ consteval size_t CompactMatrixBlock<M, Row, Col>::getSizeAtCompile() noexcept {
        return Row == 1 ? Col : Row;
    }

    template<Matrix M, size_t Row, size_t Col>
    class CompactMatrixBlock<M, Row, Col> : public LValueMatrix<CompactMatrixBlock<M, Row, Col>> {
        using This = CompactMatrixBlock<M, Row, Col>;
        using Base = LValueMatrix<This>;
        using MaybeRowCount = std::conditional_t<Row == Dynamic, size_t, Empty>;
        using MaybeColCount = std::conditional_t<Col == Dynamic, size_t, Empty>;
    private:
        decay_rvalue_t<M> mat;
        size_t fromRow;
        [[no_unique_address]] MaybeRowCount rowCount;
        size_t fromCol;
        [[no_unique_address]] MaybeColCount colCount;
    public:
        CompactMatrixBlock(M&& mat_, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount);
        CompactMatrixBlock(const This&) = default;
        CompactMatrixBlock(This&&) noexcept = default;
        ~CompactMatrixBlock() = default;
        /* Operators */
        This& operator=(const This& m);
        This& operator=(This&& m) noexcept;
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize(size_t row, size_t col);

        [[nodiscard]] auto row(this auto&&, size_t r) noexcept;
        [[nodiscard]] auto col(this auto&&, size_t c) noexcept;
        [[nodiscard]] auto rows(this auto&&, size_t fromRow, size_t rowCount) noexcept;
        [[nodiscard]] auto topRows(this auto&&, size_t to) noexcept;
        [[nodiscard]] auto bottomRows(this auto&&, size_t from) noexcept;
        [[nodiscard]] auto cols(this auto&&, size_t fromCol, size_t colCount) noexcept;
        [[nodiscard]] auto leftCols(this auto&&, size_t to) noexcept;
        [[nodiscard]] auto rightCols(this auto&&, size_t from) noexcept;
        [[nodiscard]] auto topLeftCorner(this auto&&, size_t toRow, size_t toCol) noexcept;
        [[nodiscard]] auto topLeftCorner(this auto&&, size_t to) noexcept;
        [[nodiscard]] auto topRightCorner(this auto&&, size_t toRow, size_t fromCol) noexcept;
        [[nodiscard]] auto bottomLeftCorner(this auto&&, size_t fromRow, size_t toCol) noexcept;
        [[nodiscard]] auto bottomRightCorner(this auto&&, size_t fromRow, size_t fromCol) noexcept;
        [[nodiscard]] auto bottomRightCorner(this auto&&, size_t from) noexcept;
        [[nodiscard]] auto block(this auto&&, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept;

        [[nodiscard]] auto values(this auto&&) noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] auto grads(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept;
        [[nodiscard]] size_t getCol() const noexcept;
        [[nodiscard]] size_t getOrder() const noexcept;
        [[nodiscard]] auto data_ptr(this auto&& self, size_t row, size_t col) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getRowAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getColAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept;
    };

    template<Matrix M, size_t Row, size_t Col>
    CompactMatrixBlock<M, Row, Col>::CompactMatrixBlock(M&& mat_, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount)
            : mat(std::forward<M>(mat_))
            , fromRow(fromRow)
            , rowCount(rowCount)
            , fromCol(fromCol)
            , colCount(colCount) {
        Base::checkBlock(mat, fromRow, rowCount, fromCol, colCount);
        assert(Row == Dynamic || Row == rowCount);
        assert(Col == Dynamic || Col == colCount);
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
    void CompactMatrixBlock<M, Row, Col>::resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) {
        assert(row == getRow());
        assert(col == getCol());
    }

    template<Matrix M, size_t Row, size_t Col>
    auto CompactMatrixBlock<M, Row, Col>::row(this auto&& self, size_t r) noexcept {
        assert(r < self.getRow());
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        if constexpr (MatrixMajor::isRowMatrix<M>())
            return CompactMatrixBlock<M1, 1, Col>(std::forward<M1>(m), self.fromRow + r, 1, self.fromCol, self.getCol());
        else
            return LMatrixBlock<M1, 1, Col>(std::forward<M1>(m), self.fromRow + r, 1, self.fromCol, self.getCol());
    }

    template<Matrix M, size_t Row, size_t Col>
    auto CompactMatrixBlock<M, Row, Col>::col(this auto&& self, size_t c) noexcept {
        assert(c < self.getCol());
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        if constexpr (MatrixMajor::isColMatrix<M>())
            return CompactMatrixBlock<M1, Row, 1>(std::forward<M1>(m), self.fromRow, self.getRow(), self.fromCol + c, 1);
        else
            return LMatrixBlock<M1, Row, 1>(std::forward<M1>(m), self.fromRow, self.getRow(), self.fromCol + c, 1);
    }

    template<Matrix M, size_t Row, size_t Col>
    auto CompactMatrixBlock<M, Row, Col>::rows(this auto&& self, size_t fromRow, size_t rowCount) noexcept {
        Base::checkBlock(self, fromRow, rowCount, 0, self.getCol());
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return CompactMatrixBlock<M1, Row, Col>(std::forward<M1>(m), self.fromRow + fromRow, rowCount, self.fromCol, self.getCol());
    }

    template<Matrix M, size_t Row, size_t Col>
    auto CompactMatrixBlock<M, Row, Col>::topRows(this auto&& self, size_t to) noexcept {
        Base::checkBlock(self, 0, to, 0, self.getCol());
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return CompactMatrixBlock<M1, Row, Col>(std::forward<M1>(m), self.fromRow, to, self.fromCol, self.getCol());
    }

    template<Matrix M, size_t Row, size_t Col>
    auto CompactMatrixBlock<M, Row, Col>::bottomRows(this auto&& self, size_t from) noexcept {
        Base::checkBlock(self, from, self.getRow() - from, 0, self.getCol());
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return CompactMatrixBlock<M1, Row, Col>(std::forward<M1>(m), self.fromRow + from, self.getRow() - from, self.fromCol, self.getCol());
    }

    template<Matrix M, size_t Row, size_t Col>
    auto CompactMatrixBlock<M, Row, Col>::cols(this auto&& self, size_t fromCol, size_t colCount) noexcept {
        Base::checkBlock(self, 0, self.getRow(), fromCol, colCount);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return CompactMatrixBlock<M1, Row, Col>(std::forward<M1>(m), self.fromRow, self.getRow(), self.fromCol + fromCol, colCount);
    }

    template<Matrix M, size_t Row, size_t Col>
    auto CompactMatrixBlock<M, Row, Col>::leftCols(this auto&& self, size_t to) noexcept {
        Base::checkBlock(self, 0, self.getRow(), 0, to);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return CompactMatrixBlock<M1, Row, Col>(std::forward<M1>(m), self.fromRow, self.getRow(), self.fromCol, to);
    }

    template<Matrix M, size_t Row, size_t Col>
    auto CompactMatrixBlock<M, Row, Col>::rightCols(this auto&& self, size_t from) noexcept {
        Base::checkBlock(self, 0, self.getRow(), from, self.getCol() - from);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return CompactMatrixBlock<M1, Row, Col>(std::forward<M1>(m), self.fromRow, self.getRow(), self.fromCol + from, self.getCol() - from);
    }

    template<Matrix M, size_t Row, size_t Col>
    auto CompactMatrixBlock<M, Row, Col>::topLeftCorner(this auto&& self, size_t toRow, size_t toCol) noexcept {
        Base::checkBlock(self, 0, toRow, 0, toCol);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return CompactMatrixBlock<M1, Row, Col>(std::forward<M1>(m), self.fromRow, toRow, self.fromCol, toCol);
    }

    template<Matrix M, size_t Row, size_t Col>
    auto CompactMatrixBlock<M, Row, Col>::topLeftCorner(this auto&& self, size_t to) noexcept {
        Base::checkBlock(self, 0, to, 0, to);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return CompactMatrixBlock<M1, Row, Col>(std::forward<M1>(m), self.fromRow, to, self.fromCol, to);
    }

    template<Matrix M, size_t Row, size_t Col>
    auto CompactMatrixBlock<M, Row, Col>::topRightCorner(this auto&& self, size_t toRow, size_t fromCol) noexcept {
        Base::checkBlock(self, 0, toRow, fromCol, self.getCol() - fromCol);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return CompactMatrixBlock<M1, Row, Col>(std::forward<M1>(m), self.fromRow, toRow, self.fromCol + fromCol, self.getCol() - fromCol);
    }

    template<Matrix M, size_t Row, size_t Col>
    auto CompactMatrixBlock<M, Row, Col>::bottomLeftCorner(this auto&& self, size_t fromRow, size_t toCol) noexcept {
        Base::checkBlock(self, fromRow, self.getRow() - fromRow, 0, toCol);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return CompactMatrixBlock<M1, Row, Col>(std::forward<M1>(m), self.fromRow + fromRow, self.getRow() - fromRow, self.fromCol, toCol);
    }

    template<Matrix M, size_t Row, size_t Col>
    auto CompactMatrixBlock<M, Row, Col>::bottomRightCorner(this auto&& self, size_t fromRow, size_t fromCol) noexcept {
        Base::checkBlock(self, fromRow, self.getRow() - fromRow, fromCol, self.getCol() - fromCol);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return CompactMatrixBlock<M1, Row, Col>(std::forward<M1>(m), self.fromRow + fromRow, self.getRow() - fromRow, self.fromCol + fromCol, self.getCol() - fromCol);
    }

    template<Matrix M, size_t Row, size_t Col>
    auto CompactMatrixBlock<M, Row, Col>::bottomRightCorner(this auto&& self, size_t from) noexcept {
        Base::checkBlock(self, from, self.getRow() - from, from, self.getCol() - from);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return CompactMatrixBlock<M1, Row, Col>(std::forward<M1>(m), self.fromRow + from, self.getRow() - from, self.fromCol + from, self.getCol() - from);
    }

    template<Matrix M, size_t Row, size_t Col>
    auto CompactMatrixBlock<M, Row, Col>::block(this auto&& self, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept {
        Base::checkBlock(self, fromRow, rowCount, fromCol, colCount);
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return CompactMatrixBlock<M1, Row, Col>(std::forward<M1>(m), self.fromRow + fromRow, rowCount, self.fromCol + fromCol, colCount);
    }

    template<Matrix M, size_t Row, size_t Col>
    auto CompactMatrixBlock<M, Row, Col>::values(this auto&& self) noexcept {
        auto&& v = propagate_rvalue_reference<decltype(self), M>(self.mat.values());
        using M1 = decltype(v);
        return CompactMatrixBlock<M1, Row, Col>(std::forward<M1>(v), self.fromRow, self.getRow(), self.fromCol, self.getCol());
    }

    template<Matrix M, size_t Row, size_t Col>
    template<int GradOrder>
    auto CompactMatrixBlock<M, Row, Col>::grads(this auto&& self) noexcept {
        auto&& g = propagate_rvalue_reference<decltype(self), M>(self.mat.template grads<GradOrder>());
        using M1 = decltype(g);
        return CompactMatrixBlock<M1, Row, Col>(std::forward<M1>(g), self.fromRow, self.getRow(), self.fromCol, self.getCol());
    }

    template<Matrix M, size_t Row, size_t Col>
    size_t CompactMatrixBlock<M, Row, Col>::getOrder() const noexcept {
        assert(Base::isSquare() && "[Error]: getOrder() assumes square matrix");
        return getRow();
    }

    template<Matrix M, size_t Row, size_t Col>
    auto CompactMatrixBlock<M, Row, Col>::data_ptr(this auto&& self, size_t row, size_t col) noexcept {
        assert(row < self.getRow());
        assert(col < self.getCol());
        return self.mat.data_ptr(row + self.fromRow, col + self.fromCol);
    }

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ consteval size_t CompactMatrixBlock<M, Row, Col>::getRowAtCompile() noexcept {
        return Row;
    }

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ consteval size_t CompactMatrixBlock<M, Row, Col>::getColAtCompile() noexcept {
        return Col;
    }

    template<Matrix M, size_t Row, size_t Col>
    size_t CompactMatrixBlock<M, Row, Col>::getRow() const noexcept {
        if constexpr (Row == Dynamic)
            return rowCount;
        else
            return Row;
    }

    template<Matrix M, size_t Row, size_t Col>
    size_t CompactMatrixBlock<M, Row, Col>::getCol() const noexcept {
        if constexpr (Col == Dynamic)
            return colCount;
        else
            return Col;
    }

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ consteval int CompactMatrixBlock<M, Row, Col>::getMajor() noexcept {
        return std::remove_cvref_t<M>::getMajor();
    }
}

namespace Physica {
    template<Matrix M, size_t Row, size_t Col>
    class Traits<CompactMatrixBlock<M, Row, Col>> {
    public:
        using ScalarType = std::remove_cvref_t<M>::ScalarType;
    };
}
