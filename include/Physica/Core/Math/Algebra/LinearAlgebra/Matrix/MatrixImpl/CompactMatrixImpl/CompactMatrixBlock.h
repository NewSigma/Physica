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

namespace Physica {
    template<class MatrixType, size_t Row = Dynamic, size_t Col = Dynamic> class CompactMatrixBlock;

    template<Matrix M, size_t Col>
    class CompactMatrixBlock<M, 1, Col> : public CompactVector<CompactMatrixBlock<M, 1, Col>> {
        static_assert(MatrixMajor::isRowMatrix<M>(), "[Error]: Col major does not have continuous row");
        using This = CompactMatrixBlock<M, 1, Col>;
        using Base = CompactVector<This>;
    private:
        decay_rvalue_t<M> mat;
        size_t fromRow;
        size_t fromCol;
        size_t colCount;
    public:
        CompactMatrixBlock(M&& mat_, size_t fromRow, size_t rowCount_, size_t fromCol, size_t colCount);
        CompactMatrixBlock(M&& mat_, size_t fromRow, size_t fromCol, size_t colCount);
        CompactMatrixBlock(const This&) = default;
        CompactMatrixBlock(This&&) noexcept = default;
        ~CompactMatrixBlock() = default;
        /* Operators */
        This& operator=(const This& m) { Base::operator=(m); return *this; }
        This& operator=(This&& m) noexcept { return *this = m; }
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }

        [[nodiscard]] auto&& row(this auto&&, size_t r) noexcept;

        [[nodiscard]] auto values(this auto&&) noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] auto grads(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return Col == Dynamic ? colCount : Col; }
        [[nodiscard]] auto data(this auto&& self) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return Col; }
    };

    template<Matrix M, size_t Col>
    CompactMatrixBlock<M, 1, Col>::CompactMatrixBlock(M&& mat_, size_t fromRow, [[maybe_unused]] size_t rowCount, size_t fromCol, size_t colCount)
            : CompactMatrixBlock(std::forward<M>(mat_), fromRow, fromCol, colCount) {
        assert(rowCount == 1);
    }

    template<Matrix M, size_t Col>
    CompactMatrixBlock<M, 1, Col>::CompactMatrixBlock(M&& mat_, size_t fromRow, size_t fromCol, size_t colCount)
            : mat(std::forward<M>(mat_))
            , fromRow(fromRow)
            , fromCol(fromCol)
            , colCount(colCount) {
        assert(fromRow < mat.getRow());
        assert(fromCol + colCount <= mat.getCol());
    }

    template<Matrix M, size_t Col>
    auto&& CompactMatrixBlock<M, 1, Col>::row(this auto&& self, [[maybe_unused]] size_t r) noexcept {
        assert(r == 0);
        return self;
    }

    template<Matrix M, size_t Col>
    auto CompactMatrixBlock<M, 1, Col>::values(this auto&& self) noexcept {
        auto&& v = propagate_rvalue_reference<decltype(self), M>(self.mat).values();
        using M1 = decltype(v);
        return CompactMatrixBlock<M1, 1, Col>(std::forward<M1>(v), self.fromRow, self.fromCol, self.getLength());
    }

    template<Matrix M, size_t Col>
    template<int GradOrder>
    auto CompactMatrixBlock<M, 1, Col>::grads(this auto&& self) noexcept {
        auto&& g = propagate_rvalue_reference<decltype(self), M>(self.mat).template grads<GradOrder>();
        using M1 = decltype(g);
        return CompactMatrixBlock<M1, 1, Col>(std::forward<M1>(g), self.fromRow, self.fromCol, self.getLength());
    }

    template<Matrix M, size_t Col>
    auto CompactMatrixBlock<M, 1, Col>::data(this auto&& self) noexcept {
        return self.mat.data_ptr(self.fromRow, self.fromCol);
    }

    template<Matrix M, size_t Row>
    class CompactMatrixBlock<M, Row, 1> : public CompactVector<CompactMatrixBlock<M, Row, 1>> {
        static_assert(MatrixMajor::isColMatrix<M>(), "[Error]: Row major does not have continuous col");
        using This = CompactMatrixBlock<M, Row, 1>;
        using Base = CompactVector<This>;
    private:
        decay_rvalue_t<M> mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
    public:
        CompactMatrixBlock(M&& mat_, size_t fromRow, size_t rowCount, size_t fromCol, [[maybe_unused]] size_t colCount);
        CompactMatrixBlock(M&& mat_, size_t fromRow, size_t rowCount, size_t fromCol);
        CompactMatrixBlock(const This&) = default;
        CompactMatrixBlock(This&&) noexcept = default;
        ~CompactMatrixBlock() = default;
        /* Operators */
        This& operator=(const This& m) { Base::operator=(m); return *this; }
        This& operator=(This&& m) noexcept { return *this = m; }
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }

        [[nodiscard]] auto&& col(this auto&&, size_t c) noexcept;

        [[nodiscard]] auto values(this auto&&) noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] auto grads(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return Row == Dynamic ? rowCount : Row; }
        [[nodiscard]] auto data(this auto&& self) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return Row; }
    };

    template<Matrix M, size_t Row>
    CompactMatrixBlock<M, Row, 1>::CompactMatrixBlock(M&& mat_, size_t fromRow, size_t rowCount, size_t fromCol, [[maybe_unused]] size_t colCount)
            : CompactMatrixBlock(std::forward<M>(mat_), fromRow, rowCount, fromCol) {
        assert(colCount == 1);
    }

    template<Matrix M, size_t Row>
    CompactMatrixBlock<M, Row, 1>::CompactMatrixBlock(M&& mat_, size_t fromRow, size_t rowCount, size_t fromCol)
            : mat(std::forward<M>(mat_))
            , fromRow(fromRow)
            , rowCount(rowCount)
            , fromCol(fromCol) {
        assert(fromRow + rowCount <= mat.getRow());
        assert(fromCol < mat.getCol());
    }

    template<Matrix M, size_t Row>
    auto&& CompactMatrixBlock<M, Row, 1>::col(this auto&& self, [[maybe_unused]] size_t c) noexcept {
        assert(c == 0);
        return self;
    }

    template<Matrix M, size_t Row>
    auto CompactMatrixBlock<M, Row, 1>::values(this auto&& self) noexcept {
        auto&& v = propagate_rvalue_reference<decltype(self), M>(self.mat).values();
        using M1 = decltype(v);
        return CompactMatrixBlock<M1, Row, 1>(std::forward<M1>(v), self.fromRow, self.getLength(), self.fromCol);
    }

    template<Matrix M, size_t Row>
    template<int GradOrder>
    auto CompactMatrixBlock<M, Row, 1>::grads(this auto&& self) noexcept {
        auto&& g = propagate_rvalue_reference<decltype(self), M>(self.mat).template grads<GradOrder>();
        using M1 = decltype(g);
        return CompactMatrixBlock<M1, Row, 1>(std::forward<M1>(g), self.fromRow, self.getLength(), self.fromCol);
    }

    template<Matrix M, size_t Row>
    auto CompactMatrixBlock<M, Row, 1>::data(this auto&& self) noexcept {
        return self.mat.data_ptr(self.fromRow, self.fromCol);
    }

    template<Matrix M>
    class CompactMatrixBlock<M, 1, 1> : public CompactVector<CompactMatrixBlock<M, 1, 1>> {
        using This = CompactMatrixBlock<M, 1, 1>;
        using Base = CompactVector<This>;
    private:
        decay_rvalue_t<M> mat;
        size_t fromRow;
        size_t fromCol;
    public:
        CompactMatrixBlock(M&& mat_, size_t fromRow, size_t fromCol);
        CompactMatrixBlock(const This&) = default;
        CompactMatrixBlock(This&&) noexcept = default;
        ~CompactMatrixBlock() = default;
        /* Operators */
        This& operator=(const This& m) { Base::operator=(m); return *this; }
        This& operator=(This&& m) noexcept { return *this = m; }
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t length) { assert(length == 1); }

        [[nodiscard]] auto&& row(this auto&&, size_t r) noexcept;
        [[nodiscard]] auto&& col(this auto&&, size_t c) noexcept;

        [[nodiscard]] auto values(this auto&&) noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] auto grads(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] constexpr static size_t getLength() noexcept { return 1; }
        [[nodiscard]] auto data(this auto&& self) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return 1; }
    };

    template<Matrix M>
    CompactMatrixBlock<M, 1, 1>::CompactMatrixBlock(M&& mat_, size_t fromRow, size_t fromCol)
            : mat(std::forward<M>(mat_))
            , fromRow(fromRow)
            , fromCol(fromCol) {
        assert(fromRow < mat.getRow());
        assert(fromCol < mat.getCol());
    }

    template<Matrix M>
    auto&& CompactMatrixBlock<M, 1, 1>::row(this auto&& self, [[maybe_unused]] size_t r) noexcept {
        assert(r == 0);
        return self;
    }

    template<Matrix M>
    auto&& CompactMatrixBlock<M, 1, 1>::col(this auto&& self, [[maybe_unused]] size_t c) noexcept {
        assert(c == 0);
        return self;
    }

    template<Matrix M>
    auto CompactMatrixBlock<M, 1, 1>::values(this auto&& self) noexcept {
        auto&& v = propagate_rvalue_reference<decltype(self), M>(self.mat).values();
        using M1 = decltype(v);
        return CompactMatrixBlock<M1, 1, 1>(std::forward<M1>(v), self.fromRow, self.fromCol);
    }

    template<Matrix M>
    template<int GradOrder>
    auto CompactMatrixBlock<M, 1, 1>::grads(this auto&& self) noexcept {
        auto&& g = propagate_rvalue_reference<decltype(self), M>(self.mat).template grads<GradOrder>();
        using M1 = decltype(g);
        return CompactMatrixBlock<M1, 1, 1>(std::forward<M1>(g), self.fromRow, self.fromCol);
    }

    template<Matrix M>
    [[nodiscard]] auto CompactMatrixBlock<M, 1, 1>::data(this auto&& self) noexcept {
        return self.mat.data_ptr(self.fromRow, self.fromCol);
    }

    template<Matrix M, size_t Row, size_t Col>
    class CompactMatrixBlock<M, Row, Col> : public LValueMatrix<CompactMatrixBlock<M, Row, Col>> {
        using This = CompactMatrixBlock<M, Row, Col>;
        using Base = LValueMatrix<This>;
    private:
        decay_rvalue_t<M> mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
        size_t colCount;
    public:
        CompactMatrixBlock(M&& mat_, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount);
        CompactMatrixBlock(const This&) = default;
        CompactMatrixBlock(This&&) noexcept = default;
        ~CompactMatrixBlock() = default;
        /* Operators */
        This& operator=(const This& m) { Base::operator=(m); return *this; }
        This& operator=(This&& m) noexcept { return *this = m; }
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == rowCount && col == colCount); }

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
        [[nodiscard]] size_t getRow() const noexcept { return Row == Dynamic ? rowCount : Row; }
        [[nodiscard]] size_t getCol() const noexcept { return Col == Dynamic ? colCount : Col; }
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
    auto CompactMatrixBlock<M, Row, Col>::row(this auto&& self, size_t r) noexcept {
        assert(r < self.getRow());
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        if constexpr (MatrixMajor::isRowMatrix<M>())
            return CompactMatrixBlock<M1, 1, Col>(std::forward<M1>(m), self.fromRow + r, self.fromCol, self.getCol());
        else
            return LMatrixBlock<M1, 1, Col>(std::forward<M1>(m), self.fromRow + r, self.fromCol, self.getCol());
    }

    template<Matrix M, size_t Row, size_t Col>
    auto CompactMatrixBlock<M, Row, Col>::col(this auto&& self, size_t c) noexcept {
        assert(c < self.getCol());
        auto&& m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        if constexpr (MatrixMajor::isColMatrix<M>())
            return CompactMatrixBlock<M1, Row, 1>(std::forward<M1>(m), self.fromRow, self.getRow(), self.fromCol + c);
        else
            return LMatrixBlock<M1, Row, 1>(std::forward<M1>(m), self.fromRow, self.getRow(), self.fromCol + c);
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
