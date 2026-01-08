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
        using This = ContinuousMatrixBlock<M, 1, Col>;
        using Base = ContinuousVector<This>;
    private:
        M& mat;
        size_t fromRow;
        size_t fromCol;
        size_t colCount;
    public:
        ContinuousMatrixBlock(M& mat, size_t fromRow, size_t rowCount_, size_t fromCol, size_t colCount);
        ContinuousMatrixBlock(M& mat, size_t fromRow, size_t fromCol, size_t colCount);
        ContinuousMatrixBlock(const This&) = default;
        ContinuousMatrixBlock(This&&) noexcept = default;
        ~ContinuousMatrixBlock() = default;
        /* Operators */
        This& operator=(const This& m) { Base::operator=(m); return *this; }
        This& operator=(This&& m) { return *this = m; }
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }

        [[nodiscard]] This& row([[maybe_unused]] size_t r) {
            assert(r == 0);
            return *this;
        }
        [[nodiscard]] const This& row([[maybe_unused]] size_t r) const {
            assert(r == 0);
            return *this;
        }

        [[nodiscard]] auto values(this auto&&) noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] auto grads(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return Col == Dynamic ? colCount : Col; }
        [[nodiscard]] auto data(this auto&& self) noexcept;
    };

    template<Matrix M, size_t Col>
    ContinuousMatrixBlock<M, 1, Col>::ContinuousMatrixBlock(M& mat, size_t fromRow, [[maybe_unused]] size_t rowCount, size_t fromCol, size_t colCount)
            : ContinuousMatrixBlock(mat, fromRow, fromCol, colCount) {
        assert(rowCount == 1);
    }

    template<Matrix M, size_t Col>
    ContinuousMatrixBlock<M, 1, Col>::ContinuousMatrixBlock(M& mat, size_t fromRow, size_t fromCol, size_t colCount)
            : mat(mat)
            , fromRow(fromRow)
            , fromCol(fromCol)
            , colCount(colCount) {
        assert(fromRow < mat.getRow());
        assert(fromCol + colCount <= mat.getCol());
    }

    template<Matrix M, size_t Col>
    auto ContinuousMatrixBlock<M, 1, Col>::values(this auto&& self) noexcept {
        auto&& v = self.mat.values();
        using M1 = std::remove_reference<decltype(v)>::type;
        return ContinuousMatrixBlock<M1, 1, Col>(v, self.fromRow, self.fromCol, self.colCount);
    }

    template<Matrix M, size_t Col>
    template<int GradOrder>
    auto ContinuousMatrixBlock<M, 1, Col>::grads(this auto&& self) noexcept {
        auto&& g = self.mat.grads();
        using M1 = std::remove_reference<decltype(g)>::type;
        return ContinuousMatrixBlock<M1, 1, Col>(g, self.fromRow, self.fromCol, self.colCount);
    }

    template<Matrix M, size_t Col>
    auto ContinuousMatrixBlock<M, 1, Col>::data(this auto&& self) noexcept {
        return self.mat.data_ptr(self.fromRow, self.fromCol);
    }

    template<Matrix M, size_t Row>
    class ContinuousMatrixBlock<M, Row, 1> : public ContinuousVector<ContinuousMatrixBlock<M, Row, 1>> {
        using This = ContinuousMatrixBlock<M, Row, 1>;
        using Base = ContinuousVector<This>;
    private:
        M& mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
    public:
        ContinuousMatrixBlock(M& mat, size_t fromRow, size_t rowCount, size_t fromCol, [[maybe_unused]] size_t colCount);
        ContinuousMatrixBlock(M& mat, size_t fromRow, size_t rowCount, size_t fromCol);
        ContinuousMatrixBlock(const This&) = default;
        ContinuousMatrixBlock(This&&) noexcept = default;
        ~ContinuousMatrixBlock() = default;
        /* Operators */
        This& operator=(const This& m) { Base::operator=(m); return *this; }
        This& operator=(This&& m) { return *this = m; }
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }

        [[nodiscard]] This& col([[maybe_unused]] size_t c) {
            assert(c == 0);
            return *this;
        }
        [[nodiscard]] const This& col([[maybe_unused]] size_t c) const {
            assert(c == 0);
            return *this;
        }

        [[nodiscard]] auto values(this auto&&) noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] auto grads(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return Row == Dynamic ? rowCount : Row; }
        [[nodiscard]] auto data(this auto&& self) noexcept;
    };

    template<Matrix M, size_t Row>
    ContinuousMatrixBlock<M, Row, 1>::ContinuousMatrixBlock(M& mat, size_t fromRow, size_t rowCount, size_t fromCol, [[maybe_unused]] size_t colCount)
            : ContinuousMatrixBlock(mat, fromRow, rowCount, fromCol) {
        assert(colCount == 1);
    }

    template<Matrix M, size_t Row>
    ContinuousMatrixBlock<M, Row, 1>::ContinuousMatrixBlock(M& mat, size_t fromRow, size_t rowCount, size_t fromCol)
            : mat(mat)
            , fromRow(fromRow)
            , rowCount(rowCount)
            , fromCol(fromCol) {
        assert(fromRow + rowCount <= mat.getRow());
        assert(fromCol < mat.getCol());
    }

    template<Matrix M, size_t Row>
    auto ContinuousMatrixBlock<M, Row, 1>::values(this auto&& self) noexcept {
        auto&& v = self.mat.values();
        using M1 = std::remove_reference<decltype(v)>::type;
        return ContinuousMatrixBlock<M1, Row, 1>(v, self.fromRow, self.rowCount, self.fromCol);
    }

    template<Matrix M, size_t Row>
    template<int GradOrder>
    auto ContinuousMatrixBlock<M, Row, 1>::grads(this auto&& self) noexcept {
        auto&& g = self.mat.grads();
        using M1 = std::remove_reference<decltype(g)>::type;
        return ContinuousMatrixBlock<M1, Row, 1>(g, self.fromRow, self.rowCount, self.fromCol);
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
        M& mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
    public:
        ContinuousMatrixBlock(M& mat, size_t fromRow, size_t rowCount, size_t fromCol);
        ContinuousMatrixBlock(const This&) = default;
        ContinuousMatrixBlock(This&&) noexcept = default;
        ~ContinuousMatrixBlock() = default;
        /* Operators */
        This& operator=(const This& m) { Base::operator=(m); return *this; }
        This& operator=(This&& m) { return *this = m; }
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t length) { assert(length == 1); }

        [[nodiscard]] This& row([[maybe_unused]] size_t r) {
            assert(r == 0);
            return *this;
        }
        [[nodiscard]] const This& row([[maybe_unused]] size_t r) const {
            assert(r == 0);
            return *this;
        }
        [[nodiscard]] This& col([[maybe_unused]] size_t c) {
            assert(c == 0);
            return *this;
        }
        [[nodiscard]] const This& col([[maybe_unused]] size_t c) const {
            assert(c == 0);
            return *this;
        }

        [[nodiscard]] auto values(this auto&&) noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] auto grads(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] constexpr static size_t getLength() noexcept { return 1; }
        [[nodiscard]] auto data(this auto&& self) noexcept;
    };

    template<Matrix M>
    ContinuousMatrixBlock<M, 1, 1>::ContinuousMatrixBlock(M& mat, size_t fromRow, size_t rowCount, size_t fromCol)
            : mat(mat)
            , fromRow(fromRow)
            , rowCount(rowCount)
            , fromCol(fromCol) {
        assert(fromRow + rowCount <= mat.getRow());
        assert(fromCol < mat.getCol());
    }

    template<Matrix M>
    auto ContinuousMatrixBlock<M, 1, 1>::values(this auto&& self) noexcept {
        auto&& v = self.mat.values();
        using M1 = std::remove_reference<decltype(v)>::type;
        return ContinuousMatrixBlock<M1, 1, 1>(v, self.fromRow, self.rowCount, self.fromCol);
    }

    template<Matrix M>
    template<int GradOrder>
    auto ContinuousMatrixBlock<M, 1, 1>::grads(this auto&& self) noexcept {
        auto&& g = self.mat.grads();
        using M1 = std::remove_reference<decltype(g)>::type;
        return ContinuousMatrixBlock<M1, 1, 1>(g, self.fromRow, self.rowCount, self.fromCol);
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
        M& mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
        size_t colCount;
    public:
        ContinuousMatrixBlock(M& mat, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount);
        ContinuousMatrixBlock(const This&) = default;
        ContinuousMatrixBlock(This&&) noexcept = default;
        ~ContinuousMatrixBlock() = default;
        /* Operators */
        This& operator=(const This& m) { Base::operator=(m); return *this; }
        This& operator=(This&& m) { return *this = m; }
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == rowCount && col == colCount); }

        [[nodiscard]] auto values(this auto&&) noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] auto grads(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return Row == Dynamic ? rowCount : Row; }
        [[nodiscard]] size_t getCol() const noexcept { return Col == Dynamic ? colCount : Col; }
        [[nodiscard]] auto data_ptr(this auto&& self, size_t row, size_t col) noexcept;
    };

    template<Matrix M, size_t Row, size_t Col>
    ContinuousMatrixBlock<M, Row, Col>::ContinuousMatrixBlock(M& mat, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount)
            : mat(mat)
            , fromRow(fromRow)
            , rowCount(rowCount)
            , fromCol(fromCol)
            , colCount(colCount) {
        assert(fromRow < mat.getRow());
        assert(fromCol < mat.getCol());
        assert((fromRow + rowCount) <= mat.getRow());
        assert((fromCol + colCount) <= mat.getCol());
    }

    template<Matrix M, size_t Row, size_t Col>
    auto ContinuousMatrixBlock<M, Row, Col>::values(this auto&& self) noexcept {
        auto&& v = self.mat.values();
        using M1 = std::remove_reference<decltype(v)>::type;
        return ContinuousMatrixBlock<M1, Row, Col>(v, self.fromRow, self.rowCount, self.fromCol, self.colCount);
    }

    template<Matrix M, size_t Row, size_t Col>
    template<int GradOrder>
    auto ContinuousMatrixBlock<M, Row, Col>::grads(this auto&& self) noexcept {
        auto&& g = self.mat.template grads<GradOrder>();
        using M1 = std::remove_reference<decltype(g)>::type;
        return ContinuousMatrixBlock<M1, Row, Col>(g, self.fromRow, self.rowCount, self.fromCol, self.colCount);
    }

    template<Matrix M, size_t Row, size_t Col>
    auto ContinuousMatrixBlock<M, Row, Col>::data_ptr(this auto&& self, size_t row, size_t col) noexcept {
        assert(row < self.rowCount);
        assert(col < self.colCount);
        return self.mat.data_ptr(row + self.fromRow, col + self.fromCol);
    }
}

namespace Physica {
    template<Matrix M, size_t Row, size_t Col>
    class Traits<ContinuousMatrixBlock<M, Row, Col>> {
    public:
        using ScalarType = M::ScalarType;
        constexpr static int Option = M::Option;
        constexpr static size_t RowAtCompile = Row;
        constexpr static size_t ColAtCompile = Col;
        constexpr static size_t SizeAtCompile = Row * Col;

        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = true;
    };
}
