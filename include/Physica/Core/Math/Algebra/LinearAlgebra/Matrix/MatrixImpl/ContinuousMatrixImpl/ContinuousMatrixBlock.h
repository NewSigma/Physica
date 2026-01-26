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

#include "../ContinuousMatrix.h"

namespace Physica {
    template<class MatrixType, size_t Row = Dynamic, size_t Col = Dynamic> class ContinuousMatrixBlock;

    template<Matrix M, size_t Col>
    class ContinuousMatrixBlock<M, 1, Col> : public ContinuousVector<ContinuousMatrixBlock<M, 1, Col>> {
        static_assert(Traits<M>::Option == MatrixOption::Row, "[Error]: Col major does not have continuous row");
        using This = ContinuousMatrixBlock<M, 1, Col>;
        using Base = ContinuousVector<This>;
    private:
        LazyDestroy<M> mat;
        size_t fromRow;
        size_t fromCol;
        size_t colCount;
    public:
        ContinuousMatrixBlock(M&& mat_, size_t fromRow, size_t rowCount_, size_t fromCol, size_t colCount);
        ContinuousMatrixBlock(M&& mat_, size_t fromRow, size_t fromCol, size_t colCount);
        ContinuousMatrixBlock(const This&) = default;
        ContinuousMatrixBlock(This&&) noexcept = default;
        ~ContinuousMatrixBlock() = default;
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
    };

    template<Matrix M, size_t Col>
    ContinuousMatrixBlock<M, 1, Col>::ContinuousMatrixBlock(M&& mat_, size_t fromRow, [[maybe_unused]] size_t rowCount, size_t fromCol, size_t colCount)
            : ContinuousMatrixBlock(std::forward<M>(mat_), fromRow, fromCol, colCount) {
        assert(rowCount == 1);
    }

    template<Matrix M, size_t Col>
    ContinuousMatrixBlock<M, 1, Col>::ContinuousMatrixBlock(M&& mat_, size_t fromRow, size_t fromCol, size_t colCount)
            : mat(std::forward<M>(mat_))
            , fromRow(fromRow)
            , fromCol(fromCol)
            , colCount(colCount) {
        assert(fromRow < mat.getRow());
        assert(fromCol + colCount <= mat.getCol());
    }

    template<Matrix M, size_t Col>
    auto&& ContinuousMatrixBlock<M, 1, Col>::row(this auto&& self, [[maybe_unused]] size_t r) noexcept {
        assert(r == 0);
        return self;
    }

    template<Matrix M, size_t Col>
    auto ContinuousMatrixBlock<M, 1, Col>::values(this auto&& self) noexcept {
        decltype(auto) v = propagate_rvalue_reference<decltype(self), M>(self.mat.values());
        using M1 = decltype(v);
        return ContinuousMatrixBlock<M1, 1, Col>(v, self.fromRow, self.fromCol, self.getLength());
    }

    template<Matrix M, size_t Col>
    template<int GradOrder>
    auto ContinuousMatrixBlock<M, 1, Col>::grads(this auto&& self) noexcept {
        decltype(auto) g = propagate_rvalue_reference<decltype(self), M>(self.mat.template grads<GradOrder>());
        using M1 = decltype(g);
        return ContinuousMatrixBlock<M1, 1, Col>(g, self.fromRow, self.fromCol, self.getLength());
    }

    template<Matrix M, size_t Col>
    auto ContinuousMatrixBlock<M, 1, Col>::data(this auto&& self) noexcept {
        return self.mat.data_ptr(self.fromRow, self.fromCol);
    }

    template<Matrix M, size_t Row>
    class ContinuousMatrixBlock<M, Row, 1> : public ContinuousVector<ContinuousMatrixBlock<M, Row, 1>> {
        static_assert(Traits<M>::Option == MatrixOption::Col, "[Error]: Row major does not have continuous col");
        using This = ContinuousMatrixBlock<M, Row, 1>;
        using Base = ContinuousVector<This>;
    private:
        LazyDestroy<M> mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
    public:
        ContinuousMatrixBlock(M&& mat_, size_t fromRow, size_t rowCount, size_t fromCol, [[maybe_unused]] size_t colCount);
        ContinuousMatrixBlock(M&& mat_, size_t fromRow, size_t rowCount, size_t fromCol);
        ContinuousMatrixBlock(const This&) = default;
        ContinuousMatrixBlock(This&&) noexcept = default;
        ~ContinuousMatrixBlock() = default;
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
    };

    template<Matrix M, size_t Row>
    ContinuousMatrixBlock<M, Row, 1>::ContinuousMatrixBlock(M&& mat_, size_t fromRow, size_t rowCount, size_t fromCol, [[maybe_unused]] size_t colCount)
            : ContinuousMatrixBlock(std::forward<M>(mat_), fromRow, rowCount, fromCol) {
        assert(colCount == 1);
    }

    template<Matrix M, size_t Row>
    ContinuousMatrixBlock<M, Row, 1>::ContinuousMatrixBlock(M&& mat_, size_t fromRow, size_t rowCount, size_t fromCol)
            : mat(std::forward<M>(mat_))
            , fromRow(fromRow)
            , rowCount(rowCount)
            , fromCol(fromCol) {
        assert(fromRow + rowCount <= mat.getRow());
        assert(fromCol < mat.getCol());
    }

    template<Matrix M, size_t Row>
    auto&& ContinuousMatrixBlock<M, Row, 1>::col(this auto&& self, [[maybe_unused]] size_t c) noexcept {
        assert(c == 0);
        return self;
    }

    template<Matrix M, size_t Row>
    auto ContinuousMatrixBlock<M, Row, 1>::values(this auto&& self) noexcept {
        decltype(auto) v = propagate_rvalue_reference<decltype(self), M>(self.mat.values());
        using M1 = decltype(v);
        return ContinuousMatrixBlock<M1, Row, 1>(v, self.fromRow, self.getLength(), self.fromCol);
    }

    template<Matrix M, size_t Row>
    template<int GradOrder>
    auto ContinuousMatrixBlock<M, Row, 1>::grads(this auto&& self) noexcept {
        decltype(auto) g = propagate_rvalue_reference<decltype(self), M>(self.mat.template grads<GradOrder>());
        using M1 = decltype(g);
        return ContinuousMatrixBlock<M1, Row, 1>(g, self.fromRow, self.getLength(), self.fromCol);
    }

    template<Matrix M, size_t Row>
    auto ContinuousMatrixBlock<M, Row, 1>::data(this auto&& self) noexcept {
        return self.mat.data_ptr(self.fromRow, self.fromCol);
    }

    template<Matrix M>
    class ContinuousMatrixBlock<M, 1, 1> : public ContinuousVector<ContinuousMatrixBlock<M, 1, 1>> {
        using This = ContinuousMatrixBlock<M, 1, 1>;
        using Base = ContinuousVector<This>;
    private:
        LazyDestroy<M> mat;
        size_t fromRow;
        size_t fromCol;
    public:
        ContinuousMatrixBlock(M&& mat_, size_t fromRow, size_t fromCol);
        ContinuousMatrixBlock(const This&) = default;
        ContinuousMatrixBlock(This&&) noexcept = default;
        ~ContinuousMatrixBlock() = default;
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
    };

    template<Matrix M>
    ContinuousMatrixBlock<M, 1, 1>::ContinuousMatrixBlock(M&& mat_, size_t fromRow, size_t fromCol)
            : mat(std::forward<M>(mat_))
            , fromRow(fromRow)
            , fromCol(fromCol) {
        assert(fromRow < mat.getRow());
        assert(fromCol < mat.getCol());
    }

    template<Matrix M>
    auto&& ContinuousMatrixBlock<M, 1, 1>::row(this auto&& self, [[maybe_unused]] size_t r) noexcept {
        assert(r == 0);
        return self;
    }

    template<Matrix M>
    auto&& ContinuousMatrixBlock<M, 1, 1>::col(this auto&& self, [[maybe_unused]] size_t c) noexcept {
        assert(c == 0);
        return self;
    }

    template<Matrix M>
    auto ContinuousMatrixBlock<M, 1, 1>::values(this auto&& self) noexcept {
        decltype(auto) v = propagate_rvalue_reference<decltype(self), M>(self.mat.values());
        using M1 = decltype(v);
        return ContinuousMatrixBlock<M1, 1, 1>(v, self.fromRow, 1, self.fromCol);
    }

    template<Matrix M>
    template<int GradOrder>
    auto ContinuousMatrixBlock<M, 1, 1>::grads(this auto&& self) noexcept {
        decltype(auto) g = propagate_rvalue_reference<decltype(self), M>(self.mat.template grads<GradOrder>());
        using M1 = decltype(g);
        return ContinuousMatrixBlock<M1, 1, 1>(g, self.fromRow, 1, self.fromCol);
    }

    template<Matrix M>
    [[nodiscard]] auto ContinuousMatrixBlock<M, 1, 1>::data(this auto&& self) noexcept {
        return self.mat.data_ptr(self.fromRow, self.fromCol);
    }

    template<Matrix M, size_t Row, size_t Col>
    class ContinuousMatrixBlock<M, Row, Col> : public LValueMatrix<ContinuousMatrixBlock<M, Row, Col>> {
        using This = ContinuousMatrixBlock<M, Row, Col>;
        using Base = LValueMatrix<This>;
    private:
        LazyDestroy<M> mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
        size_t colCount;
    public:
        ContinuousMatrixBlock(M&& mat_, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount);
        ContinuousMatrixBlock(const This&) = default;
        ContinuousMatrixBlock(This&&) noexcept = default;
        ~ContinuousMatrixBlock() = default;
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
    };

    template<Matrix M, size_t Row, size_t Col>
    ContinuousMatrixBlock<M, Row, Col>::ContinuousMatrixBlock(M&& mat_, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount)
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
    auto ContinuousMatrixBlock<M, Row, Col>::row(this auto&& self, size_t r) noexcept {
        assert(r < self.getRow());
        decltype(auto) m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        if constexpr (MatrixOption::isRowMatrix<M>())
            return ContinuousMatrixBlock<M1, 1, Col>(std::forward<M1>(m), self.fromRow + r, self.fromCol, self.getCol());
        else
            return LMatrixBlock<M1, 1, Col>(std::forward<M1>(m), self.fromRow + r, self.fromCol, self.getCol());
    }

    template<Matrix M, size_t Row, size_t Col>
    auto ContinuousMatrixBlock<M, Row, Col>::col(this auto&& self, size_t c) noexcept {
        assert(c < self.getCol());
        decltype(auto) m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        if constexpr (MatrixOption::isColMatrix<M>())
            return ContinuousMatrixBlock<M1, Row, 1>(std::forward<M1>(m), self.fromRow, self.getRow(), self.fromCol + c);
        else
            return LMatrixBlock<M1, Row, 1>(std::forward<M1>(m), self.fromRow, self.getRow(), self.fromCol + c);
    }

    template<Matrix M, size_t Row, size_t Col>
    auto ContinuousMatrixBlock<M, Row, Col>::rows(this auto&& self, size_t fromRow, size_t rowCount) noexcept {
        Base::checkBlock(self, fromRow, rowCount, 0, self.getCol());
        decltype(auto) m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return ContinuousMatrixBlock<M1, Row, Col>(std::forward<M1>(m), self.fromRow + fromRow, rowCount, self.fromCol, self.getCol());
    }

    template<Matrix M, size_t Row, size_t Col>
    auto ContinuousMatrixBlock<M, Row, Col>::topRows(this auto&& self, size_t to) noexcept {
        Base::checkBlock(self, 0, to, 0, self.getCol());
        decltype(auto) m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return ContinuousMatrixBlock<M1, Row, Col>(std::forward<M1>(m), self.fromRow, to, self.fromCol, self.getCol());
    }

    template<Matrix M, size_t Row, size_t Col>
    auto ContinuousMatrixBlock<M, Row, Col>::bottomRows(this auto&& self, size_t from) noexcept {
        Base::checkBlock(self, from, self.getRow() - from, 0, self.getCol());
        decltype(auto) m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return ContinuousMatrixBlock<M1, Row, Col>(std::forward<M1>(m), self.fromRow + from, self.getRow() - from, self.fromCol, self.getCol());
    }

    template<Matrix M, size_t Row, size_t Col>
    auto ContinuousMatrixBlock<M, Row, Col>::cols(this auto&& self, size_t fromCol, size_t colCount) noexcept {
        Base::checkBlock(self, 0, self.getRow(), fromCol, colCount);
        decltype(auto) m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return ContinuousMatrixBlock<M1, Row, Col>(std::forward<M1>(m), self.fromRow, self.getRow(), self.fromCol + fromCol, colCount);
    }

    template<Matrix M, size_t Row, size_t Col>
    auto ContinuousMatrixBlock<M, Row, Col>::leftCols(this auto&& self, size_t to) noexcept {
        Base::checkBlock(self, 0, self.getRow(), 0, to);
        decltype(auto) m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return ContinuousMatrixBlock<M1, Row, Col>(std::forward<M1>(m), self.fromRow, self.getRow(), self.fromCol, to);
    }

    template<Matrix M, size_t Row, size_t Col>
    auto ContinuousMatrixBlock<M, Row, Col>::rightCols(this auto&& self, size_t from) noexcept {
        Base::checkBlock(self, 0, self.getRow(), from, self.getCol() - from);
        decltype(auto) m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return ContinuousMatrixBlock<M1, Row, Col>(std::forward<M1>(m), self.fromRow, self.getRow(), self.fromCol + from, self.getCol() - from);
    }

    template<Matrix M, size_t Row, size_t Col>
    auto ContinuousMatrixBlock<M, Row, Col>::topLeftCorner(this auto&& self, size_t toRow, size_t toCol) noexcept {
        Base::checkBlock(self, 0, toRow, 0, toCol);
        decltype(auto) m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return ContinuousMatrixBlock<M1, Row, Col>(std::forward<M1>(m), self.fromRow, toRow, self.fromCol, toCol);
    }

    template<Matrix M, size_t Row, size_t Col>
    auto ContinuousMatrixBlock<M, Row, Col>::topLeftCorner(this auto&& self, size_t to) noexcept {
        Base::checkBlock(self, 0, to, 0, to);
        decltype(auto) m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return ContinuousMatrixBlock<M1, Row, Col>(std::forward<M1>(m), self.fromRow, to, self.fromCol, to);
    }

    template<Matrix M, size_t Row, size_t Col>
    auto ContinuousMatrixBlock<M, Row, Col>::topRightCorner(this auto&& self, size_t toRow, size_t fromCol) noexcept {
        Base::checkBlock(self, 0, toRow, fromCol, self.getCol() - fromCol);
        decltype(auto) m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return ContinuousMatrixBlock<M1, Row, Col>(std::forward<M1>(m), self.fromRow, toRow, self.fromCol + fromCol, self.getCol() - fromCol);
    }

    template<Matrix M, size_t Row, size_t Col>
    auto ContinuousMatrixBlock<M, Row, Col>::bottomLeftCorner(this auto&& self, size_t fromRow, size_t toCol) noexcept {
        Base::checkBlock(self, fromRow, self.getRow() - fromRow, 0, toCol);
        decltype(auto) m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return ContinuousMatrixBlock<M1, Row, Col>(std::forward<M1>(m), self.fromRow + fromRow, self.getRow() - fromRow, self.fromCol, toCol);
    }

    template<Matrix M, size_t Row, size_t Col>
    auto ContinuousMatrixBlock<M, Row, Col>::bottomRightCorner(this auto&& self, size_t fromRow, size_t fromCol) noexcept {
        Base::checkBlock(self, fromRow, self.getRow() - fromRow, fromCol, self.getCol() - fromCol);
        decltype(auto) m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return ContinuousMatrixBlock<M1, Row, Col>(std::forward<M1>(m), self.fromRow + fromRow, self.getRow() - fromRow, self.fromCol + fromCol, self.getCol() - fromCol);
    }

    template<Matrix M, size_t Row, size_t Col>
    auto ContinuousMatrixBlock<M, Row, Col>::bottomRightCorner(this auto&& self, size_t from) noexcept {
        Base::checkBlock(self, from, self.getRow() - from, from, self.getCol() - from);
        decltype(auto) m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return ContinuousMatrixBlock<M1, Row, Col>(std::forward<M1>(m), self.fromRow + from, self.getRow() - from, self.fromCol + from, self.getCol() - from);
    }

    template<Matrix M, size_t Row, size_t Col>
    auto ContinuousMatrixBlock<M, Row, Col>::block(this auto&& self, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept {
        Base::checkBlock(self, fromRow, rowCount, fromCol, colCount);
        decltype(auto) m = propagate_rvalue_reference<decltype(self), M>(self.mat);
        using M1 = decltype(m);
        return ContinuousMatrixBlock<M1, Row, Col>(std::forward<M1>(m), self.fromRow + fromRow, rowCount, self.fromCol + fromCol, colCount);
    }

    template<Matrix M, size_t Row, size_t Col>
    auto ContinuousMatrixBlock<M, Row, Col>::values(this auto&& self) noexcept {
        decltype(auto) v = propagate_rvalue_reference<decltype(self), M>(self.mat.values());
        using M1 = decltype(v);
        return ContinuousMatrixBlock<M1, Row, Col>(std::forward<M1>(v), self.fromRow, self.getRow(), self.fromCol, self.getCol());
    }

    template<Matrix M, size_t Row, size_t Col>
    template<int GradOrder>
    auto ContinuousMatrixBlock<M, Row, Col>::grads(this auto&& self) noexcept {
        decltype(auto) g = propagate_rvalue_reference<decltype(self), M>(self.mat.template grads<GradOrder>());
        using M1 = decltype(g);
        return ContinuousMatrixBlock<M1, Row, Col>(std::forward<M1>(g), self.fromRow, self.getRow(), self.fromCol, self.getCol());
    }

    template<Matrix M, size_t Row, size_t Col>
    auto ContinuousMatrixBlock<M, Row, Col>::data_ptr(this auto&& self, size_t row, size_t col) noexcept {
        assert(row < self.getRow());
        assert(col < self.getCol());
        return self.mat.data_ptr(row + self.fromRow, col + self.fromCol);
    }
}

namespace Physica {
    template<Matrix M, size_t Row, size_t Col>
    class Traits<ContinuousMatrixBlock<M, Row, Col>> {
        using M1 = std::remove_cvref<M>::type;
    public:
        using ScalarType = M1::ScalarType;
        constexpr static int Option = M1::Option;
        constexpr static size_t RowAtCompile = Row;
        constexpr static size_t ColAtCompile = Col;
        constexpr static size_t SizeAtCompile = Row * Col;

        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = true;
    };
}
